/* Author: Raghunandan Pratoori, 2021 */


#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/constrained_linear_operator.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>
#include <deal.II/physics/transformations.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <future>

namespace PhaseField
{
    using namespace dealii;

    typedef TrilinosWrappers::MPI::Vector vectorType;
    typedef TrilinosWrappers::SparseMatrix matrixType;

    
    // Input parameters from an external file
    namespace Parameters
    {
        //////////////// Class to store polynomial degree and quadrature order
        struct FESystem
        {
            // Polynomial degree and quadrature order
            unsigned int poly_degree;
            unsigned int quad_order;
            unsigned int poly_degree_c;
            unsigned int quad_order_c;

            // Function to declare all the parameters to read from input files
            static void declare_parameters(ParameterHandler &prm);

            // Function to read values of the declared parameters
            void parse_parameters(ParameterHandler &prm);
        };

        void FESystem::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Finite element system");
            {
                prm.declare_entry("Polynomial degree - displacement", "1", Patterns::Integer(0),
                                "Displacement system polynomial order");
                prm.declare_entry("Quadrature order - displacement", "2", Patterns::Integer(0),
                                "Displacement system Gauss quadrature order");
                prm.declare_entry("Polynomial degree - concentration", "1", Patterns::Integer(0),
                                "Concentration system polynomial order");
                prm.declare_entry("Quadrature order - concentration", "2", Patterns::Integer(0),
                                "Concentration system Gauss quadrature order");
            }
            prm.leave_subsection();
        }

        void FESystem::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Finite element system");
            {
                poly_degree = prm.get_integer("Polynomial degree - displacement");
                quad_order  = prm.get_integer("Quadrature order - displacement");
                poly_degree_c = prm.get_integer("Polynomial degree - concentration");
                quad_order_c  = prm.get_integer("Quadrature order - concentration");
            }
            prm.leave_subsection();
        }


        /////////////// Class to store geometry parameters
        struct Geometry
        {
            // Number of stages of refinement
            unsigned int refinement;

            static void declare_parameters(ParameterHandler &prm);

            void parse_parameters(ParameterHandler &prm);
        };

        void Geometry::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Geometry");
            {
                prm.declare_entry("Global refinement", "4", Patterns::Integer(0),
                                "Global refinement level");
            }
            prm.leave_subsection();
        }

        void Geometry::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Geometry");
            {
                refinement = prm.get_integer("Global refinement");
            }
            prm.leave_subsection();
        }


        //////////////// Class to store material parameters
        struct Materials
        {
            // Elastic properties for austenite - hexagonal system
            // 0 - 11, 1 - 12, 2 - 13, 3 - 33, 4 - 44
            std::vector<double> C_A_in;
            // Elastic properties for martensite - hexagonal system
            std::vector<double> C_M_in;
            // Lattice parameters for austenite - hexagonal system
            // 0 - a_alpha, 1 - c_alpha, 2 - a_omega, 3 - c_omega
            std::vector<double> lattice_param;
            // Kinetic coefficient
            double L;
            // Interaction energy 
            double A;
            // Thermal jump
            double delta_psi;
            // Limiting value for threshold
            double k;

            static void declare_parameters(ParameterHandler &prm);

            void parse_parameters(ParameterHandler &prm);
        };

        void Materials::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Material properties");
            {
                prm.declare_entry("C11 austenite", "0.0", Patterns::Double(0.0),
                                "C11 austenite");
                prm.declare_entry("C12 austenite", "0.0", Patterns::Double(0.0),
                                "C12 austenite");
                prm.declare_entry("C13 austenite", "0.0", Patterns::Double(0.0),
                                "C13 austenite");
                prm.declare_entry("C33 austenite", "0.0", Patterns::Double(0.0),
                                "C33 austenite");
                prm.declare_entry("C44 austenite", "0.0", Patterns::Double(0.0),
                                "C44 austenite");
                prm.declare_entry("C11 martensite", "0.0", Patterns::Double(0.0),
                                "C11 martensite");
                prm.declare_entry("C12 martensite", "0.0", Patterns::Double(0.0),
                                "C12 martensite");
                prm.declare_entry("C13 martensite", "0.0", Patterns::Double(0.0),
                                "C13 martensite");
                prm.declare_entry("C33 martensite", "0.0", Patterns::Double(0.0),
                                "C33 martensite");
                prm.declare_entry("C44 martensite", "0.0", Patterns::Double(0.0),
                                "C44 martensite");
                prm.declare_entry("a alpha", "0.0", Patterns::Double(0.0),
                                "a alpha");
                prm.declare_entry("c alpha", "0.0", Patterns::Double(0.0),
                                "c alpha");
                prm.declare_entry("a omega", "0.0", Patterns::Double(0.0),
                                "a omega");
                prm.declare_entry("c omega", "0.0", Patterns::Double(0.0),
                                "c omega");
                prm.declare_entry("kinetic coeff", "0.0", Patterns::Double(0.0),
                                "kinetic coeff");
                prm.declare_entry("interaction parameter", "0.0", Patterns::Double(0.0),
                                "interaction parameter");
                prm.declare_entry("thermal jump", "0.0", Patterns::Double(0.0),
                                "thermal jump");
                prm.declare_entry("limiting threshold", "0.0", Patterns::Double(0.0),
                                "limiting threshold");
            }
            prm.leave_subsection();
        }

        void Materials::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Material properties");
            {
                C_A_in.push_back(prm.get_double("C11 austenite"));
                C_A_in.push_back(prm.get_double("C12 austenite"));
                C_A_in.push_back(prm.get_double("C13 austenite"));
                C_A_in.push_back(prm.get_double("C33 austenite"));
                C_A_in.push_back(prm.get_double("C44 austenite"));
                C_M_in.push_back(prm.get_double("C11 martensite"));
                C_M_in.push_back(prm.get_double("C12 martensite"));
                C_M_in.push_back(prm.get_double("C13 martensite"));
                C_M_in.push_back(prm.get_double("C33 martensite"));
                C_M_in.push_back(prm.get_double("C44 martensite"));
                lattice_param.push_back(prm.get_double("a alpha"));
                lattice_param.push_back(prm.get_double("c alpha"));
                lattice_param.push_back(prm.get_double("a omega"));
                lattice_param.push_back(prm.get_double("c omega"));
                L       = prm.get_double("kinetic coeff");
                A       = prm.get_double("interaction parameter");
                delta_psi = prm.get_double("thermal jump");
                k       = prm.get_double("limiting threshold");
            }
            prm.leave_subsection();
        }


        /////////////////// Class to declare restart conditions
        struct Restart
        {
            bool restart;

            static void
            declare_parameters(ParameterHandler &prm);

            void
            parse_parameters(ParameterHandler &prm);
        };

        void Restart::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Restart");
            {
                prm.declare_entry("Restart", "false", Patterns::Bool(),
                                "restart a simulation");
            }
            prm.leave_subsection();
        }

        void Restart::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Restart");
            {
                restart = prm.get_bool("Restart");
            }
            prm.leave_subsection();
        }


        /////////// Class to store time parameters
        struct Time
        {
            double delta_t;
            double end_time;

            static void
            declare_parameters(ParameterHandler &prm);

            void
            parse_parameters(ParameterHandler &prm);
        };

        void Time::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Time");
            {
                prm.declare_entry("end time", "1e4", Patterns::Double(0.0),
                                "end time");
                prm.declare_entry("time step", "200", Patterns::Double(0.0),
                                "time step");
            }
            prm.leave_subsection();
        }

        void Time::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Time");
            {
                end_time = prm.get_double("end time");
                delta_t  = prm.get_double("time step");
            }
            prm.leave_subsection();
        }


        // Class consolidating all the above structures
        struct AllParameters :  public FESystem,
                                public Geometry,
                                public Materials,
                                public Restart,
                                public Time
        {
            void read_prm(const std::string &input_file);

            AllParameters();

            static void declare_parameters(ParameterHandler &prm);
            
            void parse_parameters(ParameterHandler &prm);
        };

        AllParameters::AllParameters()
        {}

        void AllParameters::read_prm(const std::string &input_file)
        {
            ParameterHandler prm;
            declare_parameters(prm);
            prm.parse_input(input_file);
            parse_parameters(prm);
        }

        void AllParameters::declare_parameters(ParameterHandler &prm)
        {
            FESystem::declare_parameters(prm);
            Geometry::declare_parameters(prm);
            Materials::declare_parameters(prm);
            Restart::declare_parameters(prm);
            Time::declare_parameters(prm);
        }
        
        void AllParameters::parse_parameters(ParameterHandler &prm)
        {
            FESystem::parse_parameters(prm);
            Geometry::parse_parameters(prm);
            Materials::parse_parameters(prm);
            Restart::parse_parameters(prm);
            Time::parse_parameters(prm);
        }
    } // namespace Parameters

    // Defining a global instance of input parameters
    Parameters::AllParameters parameters;

    // Defining a global variable for rotation
    const Tensor<2, 3, double> Rot_mat = Physics::Transformations::Rotations::rotation_matrix_3d(Point<3, double>(0.0, 0.0, 1.0), -1.0472);

    // A function to truncate double values to 5 decimal places
    double trunc(double number)
    {
        double scale = 1e-5;
        number = floor(number/scale)*scale;

        return number;
    }

    
    void move_file (const std::string &old_name, const std::string &new_name)
    {
        int error = system (("mv " + old_name + " " + new_name).c_str());

        // If the above call failed, e.g. because there is no command-line
        // available, try with internal functions.
        if (error != 0)
        {
            std::ifstream ifile(old_name.c_str());
            if (static_cast<bool>(ifile))
            {
                error = remove(new_name.c_str());
                AssertThrow (error == 0, ExcMessage(std::string ("Unable to remove file: "
                                                                + new_name
                                                                + ", although it seems to exist. "
                                                                + "The error code is "
                                                                + Utilities::to_string(error) + ".")));
            }

            error = rename(old_name.c_str(),new_name.c_str());
            AssertThrow (error == 0, ExcMessage(std::string ("Unable to rename files: ")
                                                +
                                                old_name + " -> " + new_name
                                                + ". The error code is "
                                                + Utilities::to_string(error) + "."));
        }
    }

    // Class to store time data
    class Time
    {
    public:
        Time(const double time_end, const double delta_t)
            :   timestep(0)
            ,   time_current(0.0)
            ,   time_end(time_end)
            ,   delta_t(delta_t)
        {}

        virtual ~Time() = default;

        double current() const
        {
            return time_current;
        }

        double end() const
        {
            return time_end;
        }
  
        double get_delta_t() const
        {
            return delta_t;
        }
        
        unsigned int get_timestep() const
        {
            return timestep;
        }
        
        void increment()
        {
            time_current += delta_t;
            ++timestep;
        }

    private:
        unsigned int timestep;
        double       time_current;
        const double time_end;
        const double delta_t;
    };


    ////////////////// Class to compute elastic modulus and stresses
    template <int dim>
    class MaterialConstitutive
    {
    public:
        MaterialConstitutive(/*const std::vector<double> C_A_in, const std::vector<double> C_M_in, const double A, const double k, const double delta_psi*/)
            :   det_F(1.0)
            ,   Fe(Tensor<2, dim>())
            ,   Fe_M2(Tensor<2, dim>())
            ,   Fe_M3(Tensor<2, dim>())
            ,   ge(Physics::Elasticity::StandardTensors<dim>::I)
            ,   ge_M2(Physics::Elasticity::StandardTensors<dim>::I)
            ,   ge_M3(Physics::Elasticity::StandardTensors<dim>::I)
            ,   Be(SymmetricTensor<2, dim>())
            ,   I1(0.0)
            ,   Ge(Physics::Elasticity::StandardTensors<dim>::I)
            ,   Ee(Tensor<2, dim>())
            ,   Ee_M2(Tensor<2, dim>())
            ,   Ee_M3(Tensor<2, dim>())
            ,   FeEe(Tensor<2, dim>())
            ,   FeEe_M2(Tensor<2, dim>())
            ,   FeEe_M3(Tensor<2, dim>())
            ,   EeEe(Tensor<2, dim>())
            ,   C_A(Vector<double> (9))
            ,   C_M(Vector<double> (9))
            ,   lambda_A(Vector<double> (3))
            ,   lambda_M(Vector<double> (3))
            ,   mu_A(Vector<double> (3))
            ,   mu_M(Vector<double> (3))
            ,   nu_A(Vector<double> (3))
            ,   nu_M(Vector<double> (3))
            ,   ki0(0.0)
        {}

        ~MaterialConstitutive(){}

        void update_material_data(const Tensor<2, dim> &F, const Tensor<2, dim> &F_e, const double &c_1, const double &c_2, const double &c_3)
        {
            Fe      = F_e;
            det_F   = determinant(F);
            ge      = Physics::Elasticity::Kinematics::b(Fe);
            Ge      = Physics::Elasticity::Kinematics::C(Fe);
            Ee      = Physics::Elasticity::Kinematics::E(Fe);
            I1      = trace(Ee);
            FeEe    = Fe*Ee;
            EeEe    = Ee*Ee;

            Fe_M2 = Rot_mat*Fe;
            ge_M2 = Physics::Elasticity::Kinematics::b(Fe_M2);
            Ee_M2 = Physics::Elasticity::Kinematics::E(Fe_M2);
            FeEe_M2 = Fe_M2*Ee_M2;

            Fe_M3 = Rot_mat*Fe_M2;
            ge_M3 = Physics::Elasticity::Kinematics::b(Fe_M3);
            Ee_M3 = Physics::Elasticity::Kinematics::E(Fe_M3);
            FeEe_M3 = Fe_M3*Ee_M3;

            c_total=c_1+c_2+c_3;
            c0=1-c_total;
            c1=c_1;
            c2=c_2;
            c3=c_3;

            
            // Populating 6x6 symmetric matrix for elastic properties
            std::vector<double> C_A(9);
            std::vector<double> C_M(9);

            C_A[0] = parameters.C_A_in[0];//C_A_11
            C_A[1] = parameters.C_A_in[0];//C_A_22
            C_A[2] = parameters.C_A_in[3];//C_A_33
            C_A[3] = parameters.C_A_in[4];//C_A_44
            C_A[4] = parameters.C_A_in[4];//C_A_55
            C_A[5] = (parameters.C_A_in[0]-parameters.C_A_in[1])/2;//C_A_66
            C_A[6] = parameters.C_A_in[1];//C_A_12
            C_A[7] = parameters.C_A_in[2];//C_A_13
            C_A[8] = parameters.C_A_in[2];//C_A_23

            C_M[0] = parameters.C_M_in[0];//C_M_11
            C_M[1] = parameters.C_M_in[0];//C_M_22
            C_M[2] = parameters.C_M_in[3];//C_M_33
            C_M[3] = parameters.C_M_in[4];//C_M_44
            C_M[4] = parameters.C_M_in[4];//C_M_55
            C_M[5] = (parameters.C_M_in[0]-parameters.C_M_in[1])/2;//C_M_66
            C_M[6] = parameters.C_M_in[1];//C_M_12
            C_M[7] = parameters.C_M_in[2];//C_M_13
            C_M[8] = parameters.C_M_in[2];//C_M_23
            
            lambda_A[0]= C_A[0]+C_A[8]+2*C_A[3]-(C_A[6]+C_A[7]+2*C_A[4]+2*C_A[5]);
            lambda_A[1]= C_A[1]+C_A[7]+2*C_A[4]-(C_A[6]+C_A[8]+2*C_A[3]+2*C_A[5]);
            lambda_A[2]= C_A[2]+C_A[6]+2*C_A[5]-(C_A[7]+C_A[8]+2*C_A[3]+2*C_A[4]);

            mu_A[0]= 0.5*(C_A[6]+C_A[7]-C_A[8]);
            mu_A[1]= 0.5*(C_A[6]+C_A[8]-C_A[7]);
            mu_A[2]= 0.5*(C_A[7]+C_A[8]-C_A[6]);

            nu_A[0]= 0.5*(C_A[4]+C_A[5]-C_A[3]);
            nu_A[1]= 0.5*(C_A[3]+C_A[5]-C_A[4]);
            nu_A[2]= 0.5*(C_A[3]+C_A[4]-C_A[5]);

            lambda_M[0]= C_M[0]+C_M[8]+2*C_M[3]-(C_M[6]+C_M[7]+2*C_M[4]+2*C_M[5]);
            lambda_M[1]= C_M[1]+C_M[7]+2*C_M[4]-(C_M[6]+C_M[8]+2*C_M[3]+2*C_M[5]);
            lambda_M[2]= C_M[2]+C_M[6]+2*C_M[5]-(C_M[7]+C_M[8]+2*C_M[3]+2*C_M[4]);

            mu_M[0]= 0.5*(C_M[6]+C_M[7]-C_M[8]);
            mu_M[1]= 0.5*(C_M[6]+C_M[8]-C_M[7]);
            mu_M[2]= 0.5*(C_M[7]+C_M[8]-C_M[6]);

            nu_M[0]= 0.5*(C_M[4]+C_M[5]-C_M[3]);
            nu_M[1]= 0.5*(C_M[3]+C_M[5]-C_M[4]);
            nu_M[2]= 0.5*(C_M[3]+C_M[4]-C_M[5]);

            // Threshold calculation
            ki0 = parameters.k - parameters.A;
        }

        // Compute Kirchoff stress
        SymmetricTensor<2, dim> get_tau() const
        {
            SymmetricTensor<2, dim> kirchhoff_stress;

            for (unsigned int n=0; n<dim; ++n)
                for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<=i; ++j)
                        kirchhoff_stress[i][j] += lambda_A[n]*c0*Ee[n][n]*Fe[i][n]*Fe[j][n]+
                                                    mu_A[n]*c0*(I1*Fe[i][n]*Fe[j][n]+Ee[n][n]*ge[i][j])+
                                                    2*nu_A[n]*c0*(Fe[i][n]*FeEe[j][n]+FeEe[i][n]*Fe[j][n])+
                                                    lambda_M[n]*c1*Ee[n][n]*Fe[i][n]*Fe[j][n]+
                                                    mu_M[n]*c1*(I1*Fe[i][n]*Fe[j][n]+Ee[n][n]*ge[i][j])+
                                                    2*nu_M[n]*c1*(Fe[i][n]*FeEe[j][n]+FeEe[i][n]*Fe[j][n])+
                                                    lambda_M[n]*c2*Ee_M2[n][n]*Fe_M2[i][n]*Fe_M2[j][n]+
                                                    mu_M[n]*c2*(I1*Fe_M2[i][n]*Fe_M2[j][n]+Ee_M2[n][n]*ge[i][j])+
                                                    2*nu_M[n]*c2*(Fe_M2[i][n]*FeEe_M2[j][n]+FeEe_M2[i][n]*Fe_M2[j][n])+
                                                    lambda_M[n]*c3*Ee[n][n]*Fe_M2[i][n]*Fe_M2[j][n]+
                                                    mu_M[n]*c3*(I1*Fe_M3[i][n]*Fe_M3[j][n]+Ee_M3[n][n]*ge[i][j])+
                                                    2*nu_M[n]*c3*(Fe_M3[i][n]*FeEe_M3[j][n]+FeEe_M3[i][n]*Fe_M3[j][n]);

            return kirchhoff_stress;
        }

        // Compute the total Jc
        SymmetricTensor<4, dim> get_Jc() const
        {
            SymmetricTensor<4,dim> elasticityTensor;

            for (unsigned int n=0; n<dim; ++n)
                for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<=i; ++j)
                        for (unsigned int k=0; k<dim; ++k)
                            for (unsigned int l=0; l<=k; ++l)
                                elasticityTensor[i][j][k][l] += lambda_A[n]*c0*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                        mu_A[n]*c0*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                        nu_A[n]*c0*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                        Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n])+
                                        lambda_M[n]*c1*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                        mu_M[n]*c1*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                        nu_M[n]*c1*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                        Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n])+
                                        lambda_M[n]*c2*Fe_M2[i][n]*Fe_M2[j][n]*Fe_M2[k][n]*Fe_M2[l][n]+
                                        mu_M[n]*c2*(Fe_M2[i][n]*Fe_M2[j][n]*ge_M2[k][l]+ ge_M2[i][j]*Fe_M2[k][n]*Fe_M2[l][n])+
                                        nu_M[n]*c2*(Fe_M2[i][n]*ge_M2[j][k]*Fe_M2[l][n]+ Fe_M2[j][n]*ge_M2[i][k]*Fe_M2[l][n]+
                                        Fe_M2[i][n]*ge_M2[j][l]*Fe_M2[k][n]+ Fe_M2[j][n]*ge_M2[i][l]*Fe_M2[k][n])+
                                        lambda_M[n]*c3*Fe_M3[i][n]*Fe_M3[j][n]*Fe_M3[k][n]*Fe_M3[l][n]+
                                        mu_M[n]*c3*(Fe_M3[i][n]*Fe_M3[j][n]*ge_M3[k][l]+ ge_M3[i][j]*Fe_M3[k][n]*Fe_M3[l][n])+
                                        nu_M[n]*c3*(Fe_M3[i][n]*ge_M3[j][k]*Fe_M3[l][n]+ Fe_M3[j][n]*ge_M3[i][k]*Fe_M3[l][n]+
                                        Fe_M3[i][n]*ge_M3[j][l]*Fe_M3[k][n]+ Fe_M3[j][n]*ge_M3[i][l]*Fe_M3[k][n]);
            
            return elasticityTensor;
        }

        // Compute Jc for Austenite
        SymmetricTensor<4, dim> get_Jc_A() const
        {
            SymmetricTensor<4,dim> elasticityTensor_A;

            for (unsigned int n=0; n<dim; ++n)
                for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int k=0; k<dim; ++k)
                            for (unsigned int l=0; l<dim; ++l)
                                elasticityTensor_A[i][j][k][l] +=
                                        lambda_A[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                        mu_A[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                        nu_A[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                        Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
            
            return elasticityTensor_A;
        }

        // Compute Jc for 1st Martensitic variant
        SymmetricTensor<4, dim> get_Jc_M1() const
        {
            SymmetricTensor<4,dim> elasticityTensor_M1;

                for (unsigned int n=0; n<dim; n++)
                    for (unsigned int i=0; i<dim; i++)
                        for (unsigned int j=0; j<=i; j++)
                            for (unsigned int k=0; k<dim; k++)
                                for (unsigned int l=0; l<=k; l++)
                                    elasticityTensor_M1[i][j][k][l] +=
                                            lambda_M[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                            mu_M[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                            nu_M[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                            Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
           
            return elasticityTensor_M1;
        }

        // Compute Jc for 2nd Martensitic variant
        SymmetricTensor<4, dim> get_Jc_M2() const
        {
            SymmetricTensor<4,dim> elasticityTensor_M2;

            for (unsigned int n=0; n<dim; ++n)
                for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<=i; ++j)
                        for (unsigned int k=0; k<dim; ++k)
                            for (unsigned int l=0; l<=k; ++l)
                            elasticityTensor_M2[i][j][k][l] +=
                                    lambda_M[n]*Fe_M2[i][n]*Fe_M2[j][n]*Fe_M2[k][n]*Fe_M2[l][n]+
                                    mu_M[n]*(Fe_M2[i][n]*Fe_M2[j][n]*ge_M2[k][l]+ ge_M2[i][j]*Fe_M2[k][n]*Fe_M2[l][n])+
                                    nu_M[n]*(Fe_M2[i][n]*ge_M2[j][k]*Fe_M2[l][n]+ Fe_M2[j][n]*ge_M2[i][k]*Fe_M2[l][n]+
                                    Fe_M2[i][n]*ge_M2[j][l]*Fe_M2[k][n]+ Fe_M2[j][n]*ge_M2[i][l]*Fe_M2[k][n]);

            return elasticityTensor_M2;
        }

        // Compute Jc for 3rd Martensitic variant
        SymmetricTensor<4, dim> get_Jc_M3() const
        {
            SymmetricTensor<4,dim> elasticityTensor_M3;

            for (unsigned int n=0; n<dim; ++n)
                for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<=i; ++j)
                        for (unsigned int k=0; k<dim; ++k)
                            for (unsigned int l=0; l<=k; ++l)
                            elasticityTensor_M3[i][j][k][l] +=
                                    lambda_M[n]*Fe_M3[i][n]*Fe_M3[j][n]*Fe_M3[k][n]*Fe_M3[l][n]+
                                    mu_M[n]*(Fe_M3[i][n]*Fe_M3[j][n]*ge_M3[k][l]+ ge_M3[i][j]*Fe_M3[k][n]*Fe_M3[l][n])+
                                    nu_M[n]*(Fe_M3[i][n]*ge_M3[j][k]*Fe_M3[l][n]+ Fe_M3[j][n]*ge_M3[i][k]*Fe_M3[l][n]+
                                    Fe_M3[i][n]*ge_M3[j][l]*Fe_M3[k][n]+ Fe_M3[j][n]*ge_M3[i][l]*Fe_M3[k][n]);
            
            return elasticityTensor_M3;
        }

        // Compute the driving force excluding the transformational work
        double get_driving_force_noStress () const
        {
            return  det_F*(parameters.delta_psi+parameters.A*(1-2*c_total));
        }

        const Tensor<2, dim> &get_Fe2() const
        {
            return Fe_M2;
        }

        double get_det_F() const
        {
            return det_F;
        }
        
        double get_threshold_ki0() const
        {
            return ki0;
        }

    protected:
        double det_F;
        Tensor<2, dim> Fe;
        Tensor<2, dim> Fe_M2;
        Tensor<2, dim> Fe_M3;
        SymmetricTensor<2, dim> ge;
        SymmetricTensor<2, dim> ge_M2;
        SymmetricTensor<2, dim> ge_M3;
        SymmetricTensor<2, dim> Be;
        double I1;
        SymmetricTensor<2, dim> Ge;
        Tensor<2, dim> Ee;
        Tensor<2, dim> Ee_M2;
        Tensor<2, dim> Ee_M3;
        Tensor<2, dim> FeEe;
        Tensor<2, dim> FeEe_M2;
        Tensor<2, dim> FeEe_M3;
        Tensor<2, dim> EeEe;
        Vector<double> C_A,C_M;
        Vector<double> lambda_A, lambda_M;
        Vector<double> mu_A, mu_M;
        Vector<double> nu_A, nu_M;
        double c_total;
        double c0,c1, c2, c3;
        double ki0;
    };


    ////////////////// Class to update quadrature point history
    template <int dim>
    class PointHistory
    {
    public:
        PointHistory()
            :   material(NULL)
            ,   F_inv(Physics::Elasticity::StandardTensors<dim>::I)
            ,   tau(SymmetricTensor<2, dim>())
            ,   Jc(SymmetricTensor<4, dim>())
        {}

        virtual ~PointHistory()
        {
            delete material;
            material = NULL;
        }

        void setup_lqp()
        {
            material = new MaterialConstitutive<dim>();

            update_values(Tensor<2, dim>(), double (), double (), double ()/*, Point<dim>()*/);
        }

        void update_values(const Tensor<2, dim> &Grad_u_n, double c1, double c2, double c3/*,
                            const Point<dim> q_point*/)
        {
            // Total deformation gradient
            F = Physics::Elasticity::Kinematics::F(Grad_u_n);

            // Transformation strain
            Tensor<2, dim> eps_t1, eps_t2, eps_t3;
            eps_t1[0][0] = (parameters.lattice_param[3] - parameters.lattice_param[0])/parameters.lattice_param[0];
            eps_t1[2][2] = (parameters.lattice_param[2] - parameters.lattice_param[1])/parameters.lattice_param[1];
            eps_t1[1][1] = (2.*sqrt(3.)*parameters.lattice_param[2] - 3.*sqrt(3.)*parameters.lattice_param[0])/(3.*sqrt(3.)*parameters.lattice_param[0]);
            eps_t2 = Rot_mat*eps_t1*invert(Rot_mat);
            eps_t3 = Rot_mat*eps_t2*invert(Rot_mat);

            // Transformation deformation gradient
            Ft = Physics::Elasticity::StandardTensors<dim>::I + eps_t1*c1 + eps_t2*c2 + eps_t3*c3;
            // Elastic deformation gradient
            Fe = F*invert(Ft);

            material -> update_material_data(F, Fe, c1, c2, c3);

            // Lagrangian strain
            E   = Physics::Elasticity::Kinematics::E(F);
            Ee  = Physics::Elasticity::Kinematics::E(Fe);
            
            F_inv = invert(F);
            F_inv_tr = transpose(F_inv);
            // Kirchoff stress
            tau = material -> get_tau();
            // Jacobian
            Jc = material -> get_Jc();
            // Driving force with no stress
            driving_force_noStress = material -> get_driving_force_noStress();
            // Threshold
            ki0 = material -> get_threshold_ki0();

            const Tensor<2, dim> temp_tensor = F_inv * Tensor<2, dim>(tau);
            const Tensor<2, dim> temp_tensor1 = temp_tensor * Fe;

            // Calculate driving force from Austenite to Martensite
            // Neglect change in elastic energy
            X10 = scalar_product(temp_tensor1, eps_t1) - driving_force_noStress - ki0;
            X20 = scalar_product(temp_tensor1, eps_t2) - driving_force_noStress - ki0;
            X30 = scalar_product(temp_tensor1, eps_t3) - driving_force_noStress - ki0;

            // Calculate driving force between Martensitic variants
            X12 = scalar_product(temp_tensor1, (eps_t1 - eps_t2));
            X13 = scalar_product(temp_tensor1, (eps_t1 - eps_t3));
            X23 = scalar_product(temp_tensor1, (eps_t2 - eps_t3));

            double c0 = 1 - c1 - c2 - c3;
            double slow = 1e-0;
            double softFactor = 1;

            // Kinetic equations
            if ((X10>0 && c1<1 && c0>0) || (X10<0 && c1>0 && c0<1))
                dc10 = slow * parameters.delta_t * parameters.L * X10;

            if ((X20>0 && c2<1 && c0>0) || (X20<0 && c2>0 && c0<1))
                dc20 = slow * parameters.delta_t * parameters.L * X20;

            if ((X30>0 && c3<1 && c0>0) || (X30<0 && c3>0 && c0<1))
                dc30 = slow * parameters.delta_t * parameters.L * X30;

            if  ((X12>0 && c1<1 && c2>0) || (X12<0 && c1>0 && c2<1))
                dc12 = slow * parameters.delta_t * parameters.L * X12;

            if  ((X13>0 && c1<1 && c3>0) || (X13<0 && c1>0 && c3<1))
                dc13 = slow * parameters.delta_t * parameters.L * X13;

            if  ((X23>0 && c2<1 && c3>0) || (X23<0 && c2>0 && c3<1))
                dc23 = slow * parameters.delta_t * parameters.L * X23;

            dc10 = (abs(dc10) < 1e-6)? 0 : dc10;
            dc20 = (abs(dc20) < 1e-6)? 0 : dc20;
            dc30 = (abs(dc30) < 1e-6)? 0 : dc30;
            dc12 = (abs(dc12) < 1e-6)? 0 : dc12;
            dc13 = (abs(dc13) < 1e-6)? 0 : dc13;
            dc23 = (abs(dc23) < 1e-6)? 0 : dc23;

            // double dci0 = dc10 + dc20 + dc30;
            // // AssertThrow(dci0 >= 0, ExcMessage("dci0 < 0"));

            // // c0 = trunc(c0);
            // // c1 = trunc(c1);
            // // c2 = trunc(c2);
            // // c3 = trunc(c3);
            // // dc10 = trunc(dc10);
            // // dc20 = trunc(dc20);
            // // dc30 = trunc(dc30);
            // // dc12 = trunc(dc12);
            // // dc13 = trunc(dc13);
            // // dc23 = trunc(dc23);


            // if(dc10>0)
            // {
            //     if (c0-dc10<0 || c1+dc10>1)
            //     dc10=std::min(1-c1,c0);
            // }
            // else if (dc10<0)
            // {
            //     if (c1-abs(dc10)<0 || c0+abs(dc10)>1)
            //     dc10=-std::min(c1,1-c0);
            // }

            // if(dc20>0)
            // {
            //     if (c0-dc20<0 || c2+dc20>1)
            //     dc20=std::min(1-c2,c0);
            // }
            // else if (dc20<0)
            // {
            //     if (c2-abs(dc20)<0 || c0+abs(dc20)>1)
            //     dc20=-std::min(c2,1-c0);
            // }

            // if(dc30>0)
            // {
            //     if (c0-dc30<0 || c3+dc30>1)
            //     dc30=std::min(1-c3,c0);
            // }
            // else if (dc30<0)
            // {
            //     if (c3-abs(dc30)<0 || c0+abs(dc30)>1)
            //     dc30=-std::min(c3,1-c0);
            // }

            // if(dc12>0)
            // {
            //     if (c2-dc12<0 || c1+dc12>1)
            //     dc12=std::min(1-c1,c2);
            // }
            // else if (dc12<0)
            // {
            //     if (c1-abs(dc12)<0 || c2+abs(dc12)>1)
            //     dc12=-std::min(c1,1-c2);
            // }

            // if(dc13>0)
            // {
            //     if (c3-dc13<0 || c1+dc13>1)
            //     dc13=std::min(1-c1,c3);
            // }
            // else if (dc13<0)
            // {
            //     if (c1-abs(dc13)<0 || c3+abs(dc13)>1)
            //     dc13=-std::min(c1,1-c3);
            // }

            // if(dc23>0)
            // {
            //     if (c3-dc23<0 || c2+dc23>1)
            //     dc23=std::min(1-c2,c3);
            // }
            // else if (dc23<0)
            // {
            //     if (c2-abs(dc23)<0 || c3+abs(dc23)>1)
            //     dc23=-std::min(c2,1-c3);
            // }

            // // Total change in concentration for each variant
            // double dc1i = dc10 + dc12 + dc13;
            // double dc2i = dc20 - dc12 + dc23;
            // double dc3i = dc30 - dc13 - dc23;

            const double dc10_old = dc10;
            const double dc20_old = dc20;
            const double dc30_old = dc30;
            const double dc12_old = dc12;
            const double dc13_old = dc13;
            const double dc23_old = dc23;


            // Constraints on change in concentration
            if ((c0 - dc10 - dc20 - dc30) < -1e-5)
            {
                dc10 = softFactor * dc10 * c0 / (abs(dc10_old) + abs(dc20_old) + abs(dc30_old));
                dc20 = softFactor * dc20 * c0 / (abs(dc10_old) + abs(dc20_old) + abs(dc30_old));
                dc30 = softFactor * dc30 * c0 / (abs(dc10_old) + abs(dc20_old) + abs(dc30_old));
            }

            if ((c1 + dc10 + dc12 + dc13) < 0)
            {
                dc10 = softFactor * dc10 * c1 / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
                dc12 = softFactor * dc12 * c1 / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
                dc13 = softFactor * dc13 * c1 / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
            }
            else if ((c1 + dc10 + dc12 + dc13) > 1)
            {
                dc10 = softFactor * dc10 * (1-c1) / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
                dc12 = softFactor * dc12 * (1-c1) / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
                dc13 = softFactor * dc13 * (1-c1) / (abs(dc10_old) + abs(dc12_old) + abs(dc13_old));
            }

            if ((c2 + dc20 - dc12 + dc23) < 0)
            {
                dc20 = softFactor * dc20 * c2 / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
                dc12 = softFactor * dc12 * c2 / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
                dc23 = softFactor * dc23 * c2 / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
            }
            else if ((c2 + dc20 - dc12 + dc23) > 1)
            {
                dc20 = softFactor * dc20 * (1-c2) / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
                dc12 = softFactor * dc12 * (1-c2) / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
                dc23 = softFactor * dc23 * (1-c2) / (abs(dc20_old) + abs(dc12_old) + abs(dc23_old));
            }

            if ((c3 + dc30 - dc13 - dc23) < 0)
            {
                dc30 = softFactor * dc30 * c3 / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
                dc13 = softFactor * dc13 * c3 / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
                dc23 = softFactor * dc23 * c3 / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
            }
            else if ((c3 + dc30 - dc13 - dc23) > 1)
            {
                dc30 = softFactor * dc30 * (1-c3) / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
                dc13 = softFactor * dc13 * (1-c3) / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
                dc23 = softFactor * dc23 * (1-c3) / (abs(dc30_old) + abs(dc13_old) + abs(dc23_old));
            }

            // // Update the change in concentration
            // dci0 = (abs(dc10 + dc20 + dc30) < 1e-6)? 0 : trunc(dc10 + dc20 + dc30);
            // dc1i = (abs(dc10 + dc12 + dc13) < 1e-6)? 0 : trunc(dc10 + dc12 + dc13);
            // dc2i = (abs(dc20 - dc12 + dc23) < 1e-6)? 0 : trunc(dc20 - dc12 + dc23);
            // dc3i = (abs(dc30 - dc13 - dc23) < 1e-6)? 0 : trunc(dc30 - dc13 - dc23);

            // double c0_up = (abs(c0 - dci0) < 1e-6)? 0 : trunc(c0 - dci0);
            // double c1_up = (abs(c1 + dc1i) < 1e-6)? 0 : trunc(c1 + dc1i);
            // double c2_up = (abs(c2 + dc2i) < 1e-6)? 0 : trunc(c2 + dc2i);
            // double c3_up = (abs(c3 + dc3i) < 1e-6)? 0 : trunc(c3 + dc3i);

            // // if (c0_up<-1e-5)
            // //     int p = 0;

            // AssertThrow(c0_up >= -1e-5, ExcMessage("c0 < 0"));
            // AssertThrow(c1_up >= -1e-5, ExcMessage("c1 < 0 in pointhistory"));
            // AssertThrow(c2_up >= -1e-5, ExcMessage("c2 < 0 in pointhistory"));
            // AssertThrow(c3_up >= -1e-5, ExcMessage("c3 < 0 in pointhistory"));
            // AssertThrow(c1_up <= 1, ExcMessage("c1 > 1 in pointhistory"));
            // AssertThrow(c2_up <= 1, ExcMessage("c2 > 1 in pointhistory"));
            // AssertThrow(c3_up <= 1, ExcMessage("c3 > 1 in pointhistory"));
        }

        const Tensor<2, dim> &get_F() const
        {
            return F;
        }

        double get_det_F() const
        {
            return material->get_det_F();
        }

        const Tensor<2, dim> &get_Fe() const
        {
            return Fe;
        }

        const Tensor<2, dim> &get_Ft() const
        {
            return Ft;
        }

        const Tensor<2, dim> &get_F_inv() const
        {
            return F_inv;
        }
        const Tensor<2, dim> &get_E() const
        {
            return E;
        }

        const Tensor<2, dim> &get_F_inv_tr() const
        {
            return F_inv_tr;
        }

        const SymmetricTensor<2, dim> &get_tau() const
        {
            return tau;
        }

        const SymmetricTensor<4, dim> &get_Jc() const
        {
            return Jc;
        }

        double get_update_c1() const
        {
            return dc10 + dc12 + dc13;
        }

        double get_update_c2() const
        {
            return dc20 - dc12 + dc23;
        }

        double get_update_c3() const
        {
            return dc30 - dc13 - dc23;
        }

    private:
        MaterialConstitutive<dim> *material;
        Tensor<2, dim> F;
        Tensor<2, dim> F_inv;
        Tensor<2, dim> F_inv_tr;
        Tensor<2, dim> Ft;
        Tensor<2, dim> E;
        SymmetricTensor<2, dim> Ee;
        SymmetricTensor<2, dim> tau;
        SymmetricTensor<4, dim> Jc;
        Tensor<2, dim> Fe;
        double driving_force_noStress;
        double ki0;
        double X10, X12, X13, X20, X23, X30;
        double dc10, dc12, dc13, dc20, dc23, dc30;
    };


    // defines the initial condition for the order parameter
    template <int dim>
    class InitialValues : public Function<dim>
    {
    public:
        InitialValues(const int &variant,const int &time_step)
            :   Function<dim>()
            ,   variant (variant)
            ,   time_step (time_step)
        {}
        virtual double value(const Point<dim>   &p,
                            const unsigned int  /*component = 0*/) const;
    private:
        const int variant;
        const int time_step;

    };

    template <int dim>
    double InitialValues<dim>::value(const Point<dim>  &p,
                                const unsigned int /*component*/) const
    {
        double tol = 1e-4;

        if (abs(p[0]-sqrt(2)*(3/8))<tol&&abs(p[1]-sqrt(2)*(3/8))<tol&&abs(p[2]-sqrt(2)*(3/8))<tol)
        {
            // return 1.0;
        }
    //   else if (time_step==0)
    return 0.0;
    }



    ////////////////// Class to run the basic FEA steps
    template <int dim>
    class Solid
    {
    public:
        Solid();

        virtual ~Solid();

        void run();
    
    private:
        void    make_grid();
        void    setup_system();
        void    determine_component_extractors();
        void    make_constraints(const int &it_nr);
        void    assemble_system();
        void    solve_nonlinear_timestep();
        unsigned int    solve();
        void    assemble_system_c();
        void    solve_c();
        void    setup_qph();
        void    update_qph_incremental();
        void    output_results() const;
        // void    output_resultant_stress();
        void    output_quad();
        void    writeQuadratureOutput(unsigned int _currentIncrement);
        void    save_checkpoint();
        void    load_triangulation();
        void    load_solution();
        void    load_time();

        MPI_Comm    mpi_communicator;
        
        parallel::distributed::Triangulation<dim>   triangulation;

        Time        time;

        ConditionalOStream               pcout;

        const FESystem<dim>         fe;
        DoFHandler<dim>             dof_handler;
        const unsigned int          dofs_per_cell;
        
        const FEValuesExtractors::Vector    u_fe;
        
        const QGauss<dim>           qf_cell;
        const QGauss<dim-1>         qf_face;
        const unsigned int          n_q_points;
        const unsigned int          n_q_points_f;
        AffineConstraints<double>   constraints;

        DoFHandler<dim>             history_dof_handler;
        FE_DGQ<dim>                 history_fe;

        std::vector<PointHistory<dim>>      quadrature_point_history;

        IndexSet                    locally_owned_dofs;
        IndexSet                    locally_relevant_dofs;

        matrixType                  tangent_matrix;
        vectorType                  system_rhs;
        vectorType                  solution;
        vectorType                  solution_update;

        FE_Q<dim>                   fe_c;
        DoFHandler<dim>             dof_handler_c;
        const unsigned int          dofs_per_cell_c;

        const QGauss<dim>           qf_cell_c;
        const unsigned int          n_q_points_c;
        AffineConstraints<double>   constraints_c;

        IndexSet                    locally_owned_dofs_c;
        IndexSet                    locally_relevant_dofs_c;

        matrixType                  mass_matrix;
        vectorType                  system_rhs_c1, system_rhs_c2, system_rhs_c3;
        vectorType                  solution_c0, solution_c1, solution_c2, solution_c3;
        vectorType                  old_solution_c1, old_solution_c2, old_solution_c3;
        vectorType                  solution_update_c1, solution_update_c2, solution_update_c3;

        // Vector<double>              resultant_cauchy_stress;
        // Vector<double>              resultant_first_piola_stress;
        // Vector<double>              resultant_second_piola_stress;
        // Vector<double>              static_cauchy_stress;
        // Vector<double>              static_first_piola_stress;
        // Vector<double>              static_second_piola_stress;
        // Vector<double>              resultant_lagrangian_strain;
        // Vector<double>              static_lagrangian_strain;
        // Vector<double>              order_parameter;
        // Vector<double>              static_order_parameter;

        std::vector<std::vector<double>>    outputQuadrature;
        // std::ofstream myfile_5;

        MappingQ<dim> mapping;

        int checkpoint_timestep;
        int checkpoint_time;
    };

    // Constructor for Solid class
    template <int dim>
    Solid<dim>::Solid()
        :   mpi_communicator(MPI_COMM_WORLD)
        ,   triangulation(mpi_communicator)
        ,   time(parameters.end_time, parameters.delta_t)
        ,   pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
        ,   fe(FE_Q<dim>(parameters.poly_degree), dim)
        ,   dof_handler(triangulation)
        ,   dofs_per_cell (fe.dofs_per_cell)
        ,   u_fe(0)
        ,   qf_cell(parameters.quad_order)
        ,   qf_face(parameters.quad_order)
        ,   n_q_points(qf_cell.size())
        ,   n_q_points_f(qf_face.size())
        ,   history_dof_handler(triangulation)
        ,   history_fe(parameters.poly_degree)
        ,   fe_c(parameters.poly_degree_c)
        ,   dof_handler_c(triangulation)
        ,   dofs_per_cell_c(fe_c.dofs_per_cell)
        ,   qf_cell_c(parameters.quad_order_c)
        ,   n_q_points_c (qf_cell_c.size())
        ,   mapping(parameters.poly_degree+1)
    {}

    // Destructor
    template <int dim>
    Solid<dim>::~Solid()
    {
        dof_handler.clear();
        dof_handler_c.clear();
    }

    // Make grid
    template <int dim>
    void Solid<dim>::make_grid()
    {
        GridGenerator::hyper_cube(triangulation, -0.5, 0.5, true);

        // Describe periodicity
        std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator>> periodicity_vector;
        GridTools::collect_periodic_faces(triangulation, 0, 1, 0, periodicity_vector);
        // GridTools::collect_periodic_faces(triangulation, 2, 3, 1, periodicity_vector);
        // GridTools::collect_periodic_faces(triangulation, 4, 5, 2, periodicity_vector);
        triangulation.add_periodicity(periodicity_vector);
        pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;

        if (parameters.restart)
        {
            load_triangulation();
        }
        else
        {
            triangulation.refine_global(1 * parameters.refinement);
        }
    }

    // Enumerate DoFs and set up matrix and vector objects to hold the data
    template <int dim>
    void Solid<dim>::setup_system()
    {
        dof_handler.distribute_dofs(fe);
        dof_handler_c.distribute_dofs(fe_c);
        history_dof_handler.distribute_dofs(history_fe);

        // DoFRenumbering::Cuthill_McKee(dof_handler);

        const unsigned int n_dofs   = dof_handler.n_dofs();
        const unsigned int n_dofs_c = dof_handler_c.n_dofs();

        // Number of active cells and DoFs
        pcout   << "   Number of active cells: "
                << triangulation.n_active_cells()
                << std::endl
                << "   Number of degrees of freedom: "
                << n_dofs + n_dofs_c
                << " (" << n_dofs << '+' << n_dofs_c << ')'
                << std::endl;

        locally_owned_dofs = dof_handler.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler,  locally_relevant_dofs);

        constraints.clear();
        constraints.reinit(locally_relevant_dofs);
        {
            DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
            // DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
            // DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, constraints);
        }
        constraints.close();
        // constraints.print(std::cout);

        DynamicSparsityPattern dsp(locally_relevant_dofs);
        DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
        Utilities::MPI::all_gather(mpi_communicator, dof_handler.locally_owned_dofs());
        SparsityTools::distribute_sparsity_pattern (dsp, dof_handler.locally_owned_dofs(),
                                                    mpi_communicator, locally_relevant_dofs);
        
        tangent_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
        solution.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        solution_update.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
        system_rhs.reinit(locally_owned_dofs, mpi_communicator);

        locally_owned_dofs_c = dof_handler_c.locally_owned_dofs();
        DoFTools::extract_locally_relevant_dofs(dof_handler_c,  locally_relevant_dofs_c);

        constraints_c.clear();
        constraints_c.reinit(locally_relevant_dofs);
        {
            // DoFTools::make_periodicity_constraints(dof_handler_c, 0, 1, 0, constraints_c);
            // DoFTools::make_periodicity_constraints(dof_handler_c, 2, 3, 1, constraints_c);
            // DoFTools::make_periodicity_constraints(dof_handler_c, 4, 5, 2, constraints_c);
        }
        constraints_c.close();

        DynamicSparsityPattern dsp_c(locally_relevant_dofs_c);
        DoFTools::make_sparsity_pattern(dof_handler_c, dsp_c, constraints_c, true);
        Utilities::MPI::all_gather(mpi_communicator, dof_handler_c.locally_owned_dofs());
        SparsityTools::distribute_sparsity_pattern (dsp_c, dof_handler_c.locally_owned_dofs(),
                                                    mpi_communicator, locally_relevant_dofs_c);
        
        mass_matrix.reinit(locally_owned_dofs_c, locally_owned_dofs_c, dsp_c, mpi_communicator);
        solution_c0.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        solution_c1.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        solution_c2.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        solution_c3.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        system_rhs_c1.reinit(locally_owned_dofs_c, mpi_communicator);
        system_rhs_c2.reinit(locally_owned_dofs_c, mpi_communicator);
        system_rhs_c3.reinit(locally_owned_dofs_c, mpi_communicator);

        if (parameters.restart)
        {
            load_solution();
        }

        // resultant_cauchy_stress.reinit(100000);
        // static_cauchy_stress.reinit   (100000);
        // resultant_first_piola_stress.reinit(100000);
        // static_first_piola_stress.reinit   (100000);
        // resultant_second_piola_stress.reinit(100000);
        // static_second_piola_stress.reinit   (100000);
        // resultant_lagrangian_strain.reinit(100000);
        // static_lagrangian_strain.reinit(100000);
        // static_order_parameter.reinit(100000);

        setup_qph();
    }
    

    // make_constraints
    template <int dim>
    void Solid<dim>::make_constraints(const int &it_nr)
    {
        if (it_nr > 1)
            return;

        constraints.clear();
        constraints.reinit (locally_relevant_dofs);

        const bool apply_dirichlet_bc = (it_nr == 0);
        const int  timestep = time.get_timestep();

        const FEValuesExtractors::Scalar x_displacement(0);
        const FEValuesExtractors::Scalar y_displacement(1);
        const FEValuesExtractors::Scalar z_displacement(2);

        // Fixing points or lines
        const double tol_boundary = 1e-4;
        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
            {
                for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)

                if     ((std::abs(cell->vertex(v)[0] - 0.) < tol_boundary) &&
                        (std::abs(cell->vertex(v)[1] - 0.) < tol_boundary) &&
                        (std::abs(cell->vertex(v)[2] - 0.) < tol_boundary))
                    {
                        constraints.add_line(cell->vertex_dof_index(v, 0));
                        constraints.add_line(cell->vertex_dof_index(v, 1));
                        constraints.add_line(cell->vertex_dof_index(v, 2));
                    }
            }
        
        {
            DoFTools::make_periodicity_constraints(dof_handler, 0, 1, 0, constraints);
            // DoFTools::make_periodicity_constraints(dof_handler, 2, 3, 1, constraints);
            // DoFTools::make_periodicity_constraints(dof_handler, 4, 5, 2, constraints);
        }
        // constraints.close();

        {
            IndexSet dofs_x0, dofs_x1;
            std::set< types::boundary_id > bid_x0 = std::set<types::boundary_id>();
            bid_x0.insert(0);
            DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(x_displacement),
                                            dofs_x0, bid_x0);
            std::set< types::boundary_id > bid_x1 = std::set<types::boundary_id>();
            bid_x1.insert(1);
            DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(x_displacement),
                                            dofs_x1, bid_x1);
            // dofs_x0.print();
            unsigned int n_dofs_x0 = dofs_x0.n_elements();
            IndexSet::ElementIterator it_x0 = dofs_x0.begin();
            IndexSet::ElementIterator it_x1 = dofs_x1.begin();

            double relative_displacement_x = 0.0;

            if(timestep<=24)
                relative_displacement_x = 5e-4;
            else if(timestep%05==0)
                relative_displacement_x = 1e-6;
            // else
            //     relative_displacement_x = 1e-4;

            for(unsigned int i = 0; i < n_dofs_x0; i++)
            {
                constraints.add_line (*it_x0);
                constraints.set_inhomogeneity(*it_x0,(apply_dirichlet_bc ? relative_displacement_x : 0.0));
                it_x0++;
                constraints.add_line (*it_x1);
                constraints.set_inhomogeneity(*it_x1,(apply_dirichlet_bc ? -2*relative_displacement_x : 0.0));
                it_x1++;
            }
        }

        // {
        //     IndexSet dofs_y0, dofs_y1;
        //     std::set< types::boundary_id > bid_y0 = std::set<types::boundary_id>();
        //     bid_y0.insert(2);
        //     DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(y_displacement),
        //                                     dofs_y0, bid_y0);
        //     std::set< types::boundary_id > bid_y1 = std::set<types::boundary_id>();
        //     bid_y1.insert(3);
        //     DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(y_displacement),
        //                                     dofs_y1, bid_y1);
        //     // dofs_y0.print();
        //     unsigned int n_dofs_y0 = dofs_y0.n_elements();
        //     IndexSet::ElementIterator it_y0 = dofs_y0.begin();
        //     IndexSet::ElementIterator it_y1 = dofs_y1.begin();

        //     double relative_displacement_y = 0.0;

        //     if(timestep<=20)
        //      relative_displacement_y = 8e-3;
        //     // else if(timestep%05==0)
        //     //  relative_displacement_y = 2e-4;
        //     else
        //      relative_displacement_y = 1e-4;

        //     for(unsigned int i = 0; i < n_dofs_y0; i++)
        //     {
        //         constraints.add_line (*it_y0);
        //         constraints.set_inhomogeneity(*it_y0,(apply_dirichlet_bc ? relative_displacement_y : 0.0));
        //         it_y0++;
        //         constraints.add_line (*it_y1);
        //         constraints.set_inhomogeneity(*it_y1,(apply_dirichlet_bc ? -2*relative_displacement_y : 0.0));
        //         it_y1++;
        //     }
        // }

        // {
        //     IndexSet dofs_z0, dofs_z1;
        //     std::set< types::boundary_id > bid_z0 = std::set<types::boundary_id>();
        //     bid_z0.insert(4);

        //     DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(z_displacement),
        //                                     dofs_z0, bid_z0);
        //     std::set< types::boundary_id > bid_z1 = std::set<types::boundary_id>();
        //     bid_z1.insert(5);

        //     DoFTools::extract_boundary_dofs(dof_handler, fe.component_mask(z_displacement),
        //                                     dofs_z1, bid_z1);
        //     // dofs_z0.print();
        //     unsigned int n_dofs_z0 = dofs_z0.n_elements();
        //     IndexSet::ElementIterator it_z0 = dofs_z0.begin();
        //     IndexSet::ElementIterator it_z1 = dofs_z1.begin();

        //     double relative_displacement_z = 0.0;

        //     if(timestep<=20)
        //      relative_displacement_z = 2e-3;
        //     // else if(timestep%05==0)
        //     //  relative_displacement_z = 2e-4;
        //     else
        //      relative_displacement_z = 1e-4;


        //     for(unsigned int i = 0; i < n_dofs_z0; i++)
        //     {
        //         constraints.add_line (*it_z0);
        //         constraints.set_inhomogeneity(*it_z0,(apply_dirichlet_bc ? relative_displacement_z : 0.0));
        //         it_z0++;
        //         constraints.add_line (*it_z1);
        //         constraints.set_inhomogeneity(*it_z1,(apply_dirichlet_bc ? -1*relative_displacement_z : 0.0));
        //         it_z1++;
        //     }
        // }
        constraints.close();
        // constraints.print(std::cout);
    }

    // Assemble system matrix and rhs
    template <int dim>
    void Solid<dim>::assemble_system ()
    {
        tangent_matrix = 0;
        system_rhs = 0;

        FEValues<dim> fe_values(fe, qf_cell,
                                update_values   | update_gradients |
                                update_quadrature_points | update_JxW_values);
        FEFaceValues<dim> fe_face_values(fe, qf_face,
                                        update_values         | update_quadrature_points  |
                                        update_normal_vectors | update_JxW_values);

        FullMatrix<double>   cell_matrix(dofs_per_cell, dofs_per_cell);
        Vector<double>       cell_rhs(dofs_per_cell);

        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        std::vector<double>                    Nx(dofs_per_cell);
        std::vector<Tensor<2, dim> >           grad_Nx(dofs_per_cell);
        std::vector<SymmetricTensor<2, dim> >  symm_grad_Nx(dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc;  ++cell)
            if (cell->is_locally_owned())
            {
                fe_values.reinit (cell);
                cell_matrix = 0;
                cell_rhs = 0;

                PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

                for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                {
                    const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();
                    const Tensor<2, dim> tau   = lqph[q_point].get_tau();
                    const SymmetricTensor<2, dim> symm_tau = lqph[q_point].get_tau();
                    const SymmetricTensor<4, dim> Jc = lqph[q_point].get_Jc();
                    const double JxW = fe_values.JxW(q_point);

                    for (unsigned int k=0; k<dofs_per_cell; ++k)
                    {
                        grad_Nx[k] = fe_values[u_fe].gradient(k, q_point)  * F_inv;
                        symm_grad_Nx[k] = symmetrize(grad_Nx[k]);
                    }


                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                    {
                        const unsigned int component_i = fe.system_to_component_index(i).first;
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                        {
                            const unsigned int component_j = fe.system_to_component_index(j).first;
                            cell_matrix(i, j) += symm_grad_Nx[i] * Jc * symm_grad_Nx[j] * JxW;
                            if (component_i == component_j)
                                cell_matrix(i, j) += grad_Nx[i][component_i] * tau * grad_Nx[j][component_j] * JxW;
                        }
                            cell_rhs(i) -= symm_grad_Nx[i] * symm_tau * JxW;
                    }
                }
                
                cell->get_dof_indices (local_dof_indices);
                constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                        local_dof_indices, tangent_matrix, system_rhs);
            }
        
        tangent_matrix.compress (VectorOperation::add);
        system_rhs.compress (VectorOperation::add);
    }


    // Assembling system matrix and RHS for the phase-field kinetic equation
    template <int dim>
    void Solid<dim>::assemble_system_c()
    {
        mass_matrix = 0;
        system_rhs_c1 = 0;
        system_rhs_c2 = 0;
        system_rhs_c3 = 0;

        FEValues<dim> fe_values_c (fe_c, qf_cell_c,
                                    update_values  | update_gradients |
                                    update_quadrature_points | update_JxW_values);

        FullMatrix<double>  cell_mass_matrix    (dofs_per_cell_c, dofs_per_cell_c);
        Vector<double>      cell_rhs_c1         (dofs_per_cell_c);
        Vector<double>      cell_rhs_c2         (dofs_per_cell_c);
        Vector<double>      cell_rhs_c3         (dofs_per_cell_c);

        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell_c);

        std::vector<double> phi(dofs_per_cell_c);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler_c.begin_active(),
        endc = dof_handler_c.end();
        for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
            {
                fe_values_c.reinit(cell);

                cell_mass_matrix    = 0;
                cell_rhs_c1         = 0;
                cell_rhs_c2         = 0;
                cell_rhs_c3         = 0;

                PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

                for (unsigned int q=0; q<n_q_points_c; ++q)
                {
                    const double dc1 = lqph[q].get_update_c1();
                    const double dc2 = lqph[q].get_update_c2();
                    const double dc3 = lqph[q].get_update_c3();

                    for (unsigned int k=0; k<dofs_per_cell_c; ++k)
                    {
                        phi[k] = fe_values_c.shape_value (k, q);
                    }

                    for (unsigned int i=0; i<dofs_per_cell_c; ++i)
                    {
                        for (unsigned int j=0; j<dofs_per_cell_c; ++j)
                        {
                            cell_mass_matrix(i,j) += phi[i] * phi[j] * fe_values_c.JxW(q);
                        }

                        cell_rhs_c1(i) +=  dc1 * phi[i] * fe_values_c.JxW (q);
                        cell_rhs_c2(i) +=  dc2 * phi[i] * fe_values_c.JxW (q);
                        cell_rhs_c3(i) +=  dc3 * phi[i] * fe_values_c.JxW (q);
                    }
                }


                cell->get_dof_indices(local_dof_indices);
                constraints_c.distribute_local_to_global(cell_mass_matrix, local_dof_indices, mass_matrix);
                constraints_c.distribute_local_to_global(cell_rhs_c1, local_dof_indices, system_rhs_c1);
                constraints_c.distribute_local_to_global(cell_rhs_c2, local_dof_indices, system_rhs_c2);
                constraints_c.distribute_local_to_global(cell_rhs_c3, local_dof_indices, system_rhs_c3);
            }

        mass_matrix.compress (VectorOperation::add);
        system_rhs_c1.compress (VectorOperation::add);
        system_rhs_c2.compress (VectorOperation::add);
        system_rhs_c3.compress (VectorOperation::add);
    }

    //setup_qph
    template <int dim>
    void Solid<dim>::setup_qph()
    {
        {
            unsigned int our_cells = 0;
            for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
                if (cell->is_locally_owned())
                    ++our_cells;
            triangulation.clear_user_data();
            {
                std::vector<PointHistory<dim> > tmp;
                tmp.swap (quadrature_point_history);
            }
            quadrature_point_history.resize (our_cells * n_q_points);

            unsigned int history_index = 0;
            for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
                if (cell->is_locally_owned())
                {
                    cell->set_user_pointer(&quadrature_point_history[history_index]);
                    history_index += n_q_points;
                }
            Assert(history_index == quadrature_point_history.size(), ExcInternalError());
        }

        for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
            if (cell->is_locally_owned())
            {
                PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

                Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
                Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

                for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                    lqph[q_point].setup_lqp();
            }
    }

    //update_qph_incremental
    template <int dim>
    void Solid<dim>::update_qph_incremental()
    {
        FEValues<dim> fe_values(fe, qf_cell,
                                update_values | update_gradients| update_quadrature_points);
        FEValues<dim> fe_values_c(fe_c, qf_cell,
                                    update_values | update_gradients| update_hessians);

        std::vector<Tensor<2, dim>> solution_grads_values(qf_cell.size());
        std::vector<double> solution_c1_values(qf_cell.size());
        std::vector<double> solution_c2_values(qf_cell.size());
        std::vector<double> solution_c3_values(qf_cell.size());

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        typename DoFHandler<dim>::active_cell_iterator
        cell_c = dof_handler_c.begin_active();
        for (; cell!=endc; ++cell, ++cell_c)
            if (cell->is_locally_owned())
            {
                PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

                Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
                Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

                Assert(solution_grads_values.size() == n_q_points, ExcInternalError());
                Assert(solution_c1_values.size() == n_q_points, ExcInternalError());

                fe_values.reinit(cell);
                fe_values_c.reinit(cell_c);

                fe_values[u_fe].get_function_gradients(solution,  solution_grads_values);
                fe_values_c.get_function_values(solution_c1,   solution_c1_values);
                fe_values_c.get_function_values(solution_c2,   solution_c2_values);
                fe_values_c.get_function_values(solution_c3,   solution_c3_values);

                for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
                {
                    lqph[q_point].update_values(solution_grads_values[q_point], solution_c1_values[q_point],
                                                solution_c2_values[q_point], solution_c3_values[q_point]/*,
                                                fe_values.quadrature_point(q_point)*/);
                }
            }
    }


    //output_results
    template <int dim>
    void Solid<dim>::output_results() const
    {
        DataOut<dim> data_out;

        // Output displacement and c
        std::vector<std::string> displacement_names;

        displacement_names.push_back ("x_displacement");
        displacement_names.push_back ("y_displacement");
        displacement_names.push_back ("z_displacement");

        data_out.add_data_vector(dof_handler, solution, displacement_names);
        data_out.add_data_vector(dof_handler_c, solution_c0, "c0");
        data_out.add_data_vector(dof_handler_c, solution_c1, "c1");
        data_out.add_data_vector(dof_handler_c, solution_c2, "c2");
        data_out.add_data_vector(dof_handler_c, solution_c3, "c3");

        /////////////////////////
        // Output norm of stress field
        Vector<double> norm_of_stress (triangulation.n_active_cells());
        {
            typename Triangulation<dim>::active_cell_iterator
            cell = triangulation.begin_active(),
            endc = triangulation.end();
            for (; cell!=endc; ++cell)
                if (cell->is_locally_owned())
                {
                SymmetricTensor<2,dim> accumulated_stress;
                for (unsigned int q=0; q<qf_cell.size(); ++q)
                    accumulated_stress +=
                    reinterpret_cast<PointHistory<dim>*>(cell->user_pointer())[q].get_tau();
                norm_of_stress(cell->active_cell_index())
                    = (accumulated_stress /
                        qf_cell.size()).norm();
                }
                else
                norm_of_stress(cell->active_cell_index()) = -1e+20;
        }
        data_out.add_data_vector (norm_of_stress, "norm_of_stress");

        ///////////////////////////////////////////////
        //Output stress componenets
        std::vector<std::vector<Vector<double>>> history_field_stress(dim, std::vector<Vector<double>>(dim)),
                                                local_history_values_at_qpoints_stress(dim, std::vector<Vector<double>>(dim)),
                                                local_history_fe_values_stress(dim, std::vector<Vector<double>>(dim));

        std::vector<std::vector< Vector<double>>> history_field_strain(dim, std::vector<Vector<double>>(dim)),
                                                local_history_values_at_qpoints_strain (dim, std::vector< Vector<double> >(dim)),
                                                local_history_fe_values_strain (dim, std::vector< Vector<double> >(dim));

        for (unsigned int i=0; i<dim; i++)
            for (unsigned int j=0; j<dim; j++)
            {
                history_field_stress[i][j].reinit(history_dof_handler.n_dofs());
                local_history_values_at_qpoints_stress[i][j].reinit(qf_cell.size());
                local_history_fe_values_stress[i][j].reinit(history_fe.dofs_per_cell);
            }

        for (unsigned int i=0; i<dim; i++)
            for (unsigned int j=0; j<dim; j++)
            {
                history_field_strain[i][j].reinit(history_dof_handler.n_dofs());
                local_history_values_at_qpoints_strain[i][j].reinit(qf_cell.size());
                local_history_fe_values_strain[i][j].reinit(history_fe.dofs_per_cell);
            }

        Vector<double> history_field_drivingforce_c1,
                        history_field_drivingforce_c2,
                        history_field_drivingforce_c3,
                        local_history_values_at_qpoints_drivingforce_c1,
                        local_history_values_at_qpoints_drivingforce_c2,
                        local_history_values_at_qpoints_drivingforce_c3,
                        local_history_fe_values_drivingforce_c1,
                        local_history_fe_values_drivingforce_c2,
                        local_history_fe_values_drivingforce_c3;

        history_field_drivingforce_c1.reinit(history_dof_handler.n_dofs());
        history_field_drivingforce_c2.reinit(history_dof_handler.n_dofs());
        history_field_drivingforce_c3.reinit(history_dof_handler.n_dofs());
        local_history_values_at_qpoints_drivingforce_c1.reinit(qf_cell.size());
        local_history_values_at_qpoints_drivingforce_c2.reinit(qf_cell.size());
        local_history_values_at_qpoints_drivingforce_c3.reinit(qf_cell.size());
        local_history_fe_values_drivingforce_c1.reinit(history_fe.dofs_per_cell);
        local_history_fe_values_drivingforce_c2.reinit(history_fe.dofs_per_cell);
        local_history_fe_values_drivingforce_c3.reinit(history_fe.dofs_per_cell);


        FullMatrix<double> qpoint_to_dof_matrix(history_fe.dofs_per_cell, qf_cell.size());
        FETools::compute_projection_from_quadrature_points_matrix(history_fe, qf_cell, qf_cell, qpoint_to_dof_matrix);

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end(),
        dg_cell = history_dof_handler.begin_active();
        for (; cell!=endc; ++cell, ++dg_cell)
            if (cell->is_locally_owned())
            {
                PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
                Assert (lqph >= &quadrature_point_history.front(), ExcInternalError());
                Assert (lqph < &quadrature_point_history.back(), ExcInternalError());
                
                for (unsigned int i=0; i<dim; i++)
                    for (unsigned int j=0; j<dim; j++)
                    {
                        for (unsigned int q=0; q<qf_cell.size(); ++q)
                        {
                            local_history_values_at_qpoints_stress[i][j](q) = (lqph[q].get_tau()[i][j])/(lqph[q].get_det_F());
                            qpoint_to_dof_matrix.vmult(local_history_fe_values_stress[i][j], local_history_values_at_qpoints_stress[i][j]);
                            dg_cell->set_dof_values(local_history_fe_values_stress[i][j], history_field_stress[i][j]);

                            local_history_values_at_qpoints_strain[i][j](q) = lqph[q].get_E()[i][j];
                            qpoint_to_dof_matrix.vmult(local_history_fe_values_strain[i][j], local_history_values_at_qpoints_strain[i][j]);
                            dg_cell->set_dof_values(local_history_fe_values_strain[i][j], history_field_strain[i][j]);

                            local_history_values_at_qpoints_drivingforce_c1(q) = lqph[q].get_update_c1();
                            qpoint_to_dof_matrix.vmult(local_history_fe_values_drivingforce_c1, local_history_values_at_qpoints_drivingforce_c1);
                            dg_cell->set_dof_values(local_history_fe_values_drivingforce_c1, history_field_drivingforce_c1);

                            local_history_values_at_qpoints_drivingforce_c2(q) = lqph[q].get_update_c2();
                            qpoint_to_dof_matrix.vmult(local_history_fe_values_drivingforce_c2, local_history_values_at_qpoints_drivingforce_c2);
                            dg_cell->set_dof_values(local_history_fe_values_drivingforce_c2, history_field_drivingforce_c2);

                            local_history_values_at_qpoints_drivingforce_c3(q) = lqph[q].get_update_c3();
                            qpoint_to_dof_matrix.vmult(local_history_fe_values_drivingforce_c3, local_history_values_at_qpoints_drivingforce_c3);
                            dg_cell->set_dof_values(local_history_fe_values_drivingforce_c3, history_field_drivingforce_c3);
                        }
                    }
            }

        std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation2(1, DataComponentInterpretation::component_is_scalar);

        data_out.add_data_vector(history_dof_handler, history_field_stress[0][0], "sigma_11", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_stress[1][1], "sigma_22", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_stress[2][2], "sigma_33", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_stress[0][1], "sigma_12", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_stress[0][2], "sigma_13", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_stress[1][2], "sigma_23", data_component_interpretation2);

        data_out.add_data_vector(history_dof_handler, history_field_strain[0][0], "E_11", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_strain[1][1], "E_22", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_strain[2][2], "E_33", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_strain[0][1], "E_12", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_strain[0][2], "E_13", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_strain[1][2], "E_23", data_component_interpretation2);

        data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c1, "dc1", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c2, "dc2", data_component_interpretation2);
        data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c3, "dc3", data_component_interpretation2);

        //////////////////////////
        // writing output files
        MappingQEulerian<dim, vectorType> q_mapping(parameters.poly_degree, dof_handler, solution);

        Vector<float> subdomain(triangulation.n_active_cells());
        for (unsigned int i=0; i<subdomain.size(); ++i)
            subdomain(i) = triangulation.locally_owned_subdomain();
        
        data_out.add_data_vector(subdomain, "subdomain");

        data_out.build_patches(q_mapping, parameters.poly_degree);

        const unsigned int cycle = time.get_timestep();

        const std::string filename = ("solution-" +
                                        Utilities::int_to_string (cycle, 4) +
                                        "." +
                                        Utilities::int_to_string
                                        (triangulation.locally_owned_subdomain(), 4));
        std::ofstream output((filename + ".vtu").c_str());
        data_out.write_vtu(output);

        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
            std::vector<std::string> filenames;
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
            {
                filenames.push_back("solution-" + Utilities::int_to_string (cycle, 4) +
                                    "." + Utilities::int_to_string (i, 4) + ".vtu");
            }
            std::ofstream master_output(("solution-" + Utilities::int_to_string (cycle, 4) + ".pvtu").c_str());
            data_out.write_pvtu_record (master_output, filenames);
        }
    }

    //output quadrature values
    template <int dim>
    void Solid<dim>::output_quad()
    {
        FEValues<dim> fe_values (fe, qf_cell,
                        update_values   | update_gradients |
                        update_quadrature_points | update_JxW_values);
        outputQuadrature.clear();

        typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());
        for (; cell != endc; ++cell)
        {
            if (cell->is_locally_owned())
            {
                fe_values.reinit(cell);
                // loop over quadrature points
                for (unsigned int q=0; q<n_q_points; ++q)
                {
                    std::vector<double> temp;

                    const Tensor<2, dim> F = lqph[q].get_F();
                    const double J = lqph[q].get_det_F();
                    const Tensor<2, dim> Fe = lqph[q].get_Fe();
                    const Tensor<2, dim> Ft = lqph[q].get_Ft();
                    const Tensor<2, dim> E = lqph[q].get_E();
                    const Tensor<2, dim> s = lqph[q].get_tau()/J;
                    const double dc1 = lqph[q].get_update_c1();
                    const double dc2 = lqph[q].get_update_c2();
                    const double dc3 = lqph[q].get_update_c3();
                    const double X10 = lqph[q].get_X10();
                    const double X12 = lqph[q].get_X12();
                    const double X13 = lqph[q].get_X13();
                    const double X20 = lqph[q].get_X20();
                    const double X30 = lqph[q].get_X30();
                    const double X21 = lqph[q].get_X21();
                    const double X31 = lqph[q].get_X31();
                    const double X23 = lqph[q].get_X23();
                    const double X32 = lqph[q].get_X32();
                    const double dc10 = lqph[q].get_dc10();
                    const double dc12 = lqph[q].get_dc12();
                    const double dc13 = lqph[q].get_dc13();
                    const double dc20 = lqph[q].get_dc20();
                    const double dc30 = lqph[q].get_dc30();
                    const double dc21 = lqph[q].get_dc21();
                    const double dc31 = lqph[q].get_dc31();
                    const double dc23 = lqph[q].get_dc23();
                    const double dc32 = lqph[q].get_dc32();


                    temp.push_back(fe_values.JxW(q));

                    temp.push_back(fe_values.get_quadrature_points()[q][0]);
                    temp.push_back(fe_values.get_quadrature_points()[q][1]);
                    temp.push_back(fe_values.get_quadrature_points()[q][2]);

                    temp.push_back(F[0][0]);
                    temp.push_back(F[0][1]);
                    temp.push_back(F[0][2]);
                    temp.push_back(F[1][0]);
                    temp.push_back(F[1][1]);
                    temp.push_back(F[1][2]);
                    temp.push_back(F[2][0]);
                    temp.push_back(F[2][1]);
                    temp.push_back(F[2][2]);

                    temp.push_back(J);

                    temp.push_back(Fe[0][0]);
                    temp.push_back(Fe[0][1]);
                    temp.push_back(Fe[0][2]);
                    temp.push_back(Fe[1][0]);
                    temp.push_back(Fe[1][1]);
                    temp.push_back(Fe[1][2]);
                    temp.push_back(Fe[2][0]);
                    temp.push_back(Fe[2][1]);
                    temp.push_back(Fe[2][2]);

                    temp.push_back(Ft[0][0]);
                    temp.push_back(Ft[0][1]);
                    temp.push_back(Ft[0][2]);
                    temp.push_back(Ft[1][0]);
                    temp.push_back(Ft[1][1]);
                    temp.push_back(Ft[1][2]);
                    temp.push_back(Ft[2][0]);
                    temp.push_back(Ft[2][1]);
                    temp.push_back(Ft[2][2]);

                    temp.push_back(E[0][0]);
                    temp.push_back(E[0][1]);
                    temp.push_back(E[0][2]);
                    temp.push_back(E[1][0]);
                    temp.push_back(E[1][1]);
                    temp.push_back(E[1][2]);
                    temp.push_back(E[2][0]);
                    temp.push_back(E[2][1]);
                    temp.push_back(E[2][2]);

                    temp.push_back(s[0][0]);
                    temp.push_back(s[0][1]);
                    temp.push_back(s[0][2]);
                    temp.push_back(s[1][0]);
                    temp.push_back(s[1][1]);
                    temp.push_back(s[1][2]);
                    temp.push_back(s[2][0]);
                    temp.push_back(s[2][1]);
                    temp.push_back(s[2][2]);

                    temp.push_back(dc1);
                    temp.push_back(dc2);
                    temp.push_back(dc3);

                    temp.push_back(X10);
                    temp.push_back(X20);
                    temp.push_back(X30);
                    temp.push_back(X12);
                    temp.push_back(X13);
                    temp.push_back(X21);
                    temp.push_back(X23);
                    temp.push_back(X31);
                    temp.push_back(X32);

                    temp.push_back(dc10);
                    temp.push_back(dc20);
                    temp.push_back(dc30);
                    temp.push_back(dc12);
                    temp.push_back(dc13);
                    temp.push_back(dc21);
                    temp.push_back(dc23);
                    temp.push_back(dc31);
                    temp.push_back(dc32);

                    outputQuadrature.push_back(temp);
                }
            }
            writeQuadratureOutput(this->time.get_timestep());
        }
    }

    template <int dim>
    void Solid<dim>::writeQuadratureOutput(unsigned int _currentIncrement)
    {
        this->pcout << "writing Quadrature data to file\n";
        {
            std::string fileName("QuadratureOutputs");
            std::string fileExtension(".csv");
            fileName += std::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
            std::ofstream file((fileName + Utilities::int_to_string(_currentIncrement,4)+fileExtension).c_str());
            char buffer[200];
            if (file.is_open())
            {
                for (std::vector<std::vector<double> >::iterator it = outputQuadrature.begin(); it != outputQuadrature.end(); ++it)
                {
                    for (std::vector<double>::iterator it2 = it->begin(); it2 != it->end(); ++it2)
                    {
                        sprintf(buffer, "%8.5e ,", *it2);
                        file << buffer;
                    }
                    file << std::endl;
                }
                file.close();
            }
            else
            {
                this->pcout << "Unable to open file for writing quadrature outputs\n";
                exit(1);
            }


            std::string fileName2("QuadratureOutputs");
            std::ofstream file2((fileName2 + Utilities::int_to_string(_currentIncrement,4)+fileExtension).c_str());
            for(unsigned int proc = 0; proc<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); proc++)
            {
                std::string fileName3("QuadratureOutputs");
                fileName3 += std::to_string(proc);
                std::ofstream file3((fileName3 + Utilities::int_to_string(_currentIncrement,4)+fileExtension).c_str(), std::ofstream::in);
                file2 << file3.rdbuf();
                //delete file from processor proc
                remove((fileName3 + Utilities::int_to_string(_currentIncrement,4)+fileExtension).c_str());
            }
            file2.close();
        }
    }

    template <int dim>
    void Solid<dim>::save_checkpoint()
    {
        this->pcout << "Saving checkpoint\n";
        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        {
            static bool previous_snapshot_exists = (parameters.restart == true);

            if (previous_snapshot_exists == true)
            {
                move_file("restart.mesh","restart.mesh.old");
                move_file("restart.mesh.info","restart.mesh.info.old");
                move_file("restart.time.info","restart.time.info.old");
            }
            previous_snapshot_exists = true;
        }

        // Save triangulation and solution vectors
        parallel::distributed::SolutionTransfer<dim, vectorType> solution_transfer(dof_handler);
        parallel::distributed::SolutionTransfer<dim, vectorType> solution_transfer_c(dof_handler_c);

        std::vector<const vectorType *> displacement_transfer;
        std::vector<const vectorType *> volFraction_transfer;

        displacement_transfer.push_back(&solution);
        volFraction_transfer.push_back(&solution_c1);
        volFraction_transfer.push_back(&solution_c2);
        volFraction_transfer.push_back(&solution_c3);

        solution_transfer.prepare_for_serialization(displacement_transfer);
        solution_transfer_c.prepare_for_serialization(volFraction_transfer);

        triangulation.save("restart.mesh");

        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
        {
            std::ofstream time_info_file;
            time_info_file.open("restart.time.info");
            time_info_file << time.get_timestep() << " (currentIncrement)\n";
            time_info_file << time.current() << " (currentTime)\n";
            time_info_file.close();
        }
    }

    template <int dim>
    void Solid<dim>::load_triangulation()
    {
        this->pcout << "Loading triangulation\n";

        triangulation.load("restart.mesh");
    }

    template <int dim>
    void Solid<dim>::load_solution()
    {
        this->pcout << "Loading solution\n";

        std::vector<vectorType *> displacement_transfer(1);
        std::vector<vectorType *> volFraction_transfer(3);

        parallel::distributed::SolutionTransfer<dim, vectorType> solution_transfer(dof_handler);
        parallel::distributed::SolutionTransfer<dim, vectorType> solution_transfer_c(dof_handler_c);

        vectorType displacement_system(system_rhs);
        vectorType c1_system(system_rhs_c1); 
        vectorType c2_system(system_rhs_c2); 
        vectorType c3_system(system_rhs_c3);

        displacement_transfer[0] = &displacement_system;
        volFraction_transfer[0] = &c1_system;
        volFraction_transfer[1] = &c2_system;
        volFraction_transfer[2] = &c3_system;

        solution_transfer.deserialize(displacement_transfer);
        solution_transfer_c.deserialize(volFraction_transfer);

        solution = displacement_system;
        solution_c1 = c1_system;
        solution_c2 = c2_system;
        solution_c3 = c3_system;
    }

    template <int dim>
    void Solid<dim>::load_time()
    {
        this->pcout << "Loading time data\n";

        std::ifstream time_info_file;
        time_info_file.open("restart.time.info");
        std::string line;
        std::getline(time_info_file, line);
        line.erase(line.end()-19,line.end());
        checkpoint_timestep = Utilities::string_to_int(line);

        std::getline(time_info_file, line);
        line.erase(line.end()-14,line.end());
        checkpoint_time = Utilities::string_to_double(line);
        time_info_file.close();
    }

    // template <int dim>
    // void Solid<dim>::output_resultant_stress()
    // {

    //     FEFaceValues<dim> fe_face_values (fe, qf_face,
    //                                     update_values         | update_quadrature_points  |
    //                                     update_normal_vectors | update_JxW_values);
    //     double resultant_force = 0.0;
    //     double resultant_pseudo_force = 0.0;
    //     double current_face_area = 0.0;
    //     double referrence_face_area = 0.0;
    //     double resultant_E = 0.0;
    //     typename DoFHandler<dim>::active_cell_iterator
    //     cell = dof_handler.begin_active(),
    //     endc = dof_handler.end();
    //     for (; cell != endc; ++cell)
    //         if (cell->is_locally_owned())
    //     {

    //         PointHistory<dim> *lqph = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    //         for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    //             if (cell->face(face)->at_boundary() && (cell->face(face)->boundary_id() == 5))
    //             {
    //                 fe_face_values.reinit (cell, face);
    //                 for (unsigned int q_point=0; q_point<n_q_points_f; ++q_point)
    //                 {
    //                     const Tensor<2, dim> F_inv_tr   = lqph[q_point].get_F_inv_tr();
    //                     const Tensor<2, dim> F_inv      = lqph[q_point].get_F_inv();
    //                     const Tensor<2, dim> tau        = lqph[q_point].get_tau();
    //                     const Tensor<2, dim> E          = lqph[q_point].get_E();
    //                     const double J = lqph[q_point].get_det_F();

    //                     const Tensor<1, dim> element_force = (tau*F_inv_tr) * fe_face_values.normal_vector(q_point);
    //                     const Tensor<1, dim> element_pseudo_force = F_inv*element_force;

    //                     const double current_element_area = J * ((F_inv_tr*fe_face_values.normal_vector(q_point)).norm());

    //                     resultant_force += element_force [2]*fe_face_values.JxW(q_point);
    //                     resultant_pseudo_force += element_pseudo_force [2]*fe_face_values.JxW(q_point);
    //                     current_face_area += current_element_area* fe_face_values.JxW(q_point);
    //                     referrence_face_area += fe_face_values.JxW(q_point);
    //                     resultant_E += E[2][2]*fe_face_values.JxW(q_point);
    //                 }
    //             }
    //     }


    //     resultant_cauchy_stress[time.get_timestep()]        = -1*Utilities::MPI::sum(resultant_force, mpi_communicator )/
    //                                                         Utilities::MPI::sum(current_face_area, mpi_communicator );
    //     resultant_first_piola_stress[time.get_timestep()]   = -1*Utilities::MPI::sum(resultant_force, mpi_communicator )/
    //                                                         Utilities::MPI::sum(referrence_face_area, mpi_communicator );
    //     resultant_second_piola_stress[time.get_timestep()]  = -1*Utilities::MPI::sum(resultant_pseudo_force, mpi_communicator )/
    //                                                         Utilities::MPI::sum(referrence_face_area, mpi_communicator );
    //     resultant_lagrangian_strain[time.get_timestep()]    = -1*Utilities::MPI::sum(resultant_E, mpi_communicator )/
    //                                                         Utilities::MPI::sum(referrence_face_area, mpi_communicator );

    //     if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    //     {

    //         std::ofstream myfile_1;
    //         std::ofstream myfile_2;
    //         std::ofstream myfile_3;
    //         std::ofstream myfile_4;
            
    //         myfile_1.open ("resultant_cauchy_stress.txt");
    //         myfile_2.open ("resultant_first_piola_stress.txt");
    //         myfile_3.open ("resultant_second_piola_stress.txt");
    //         myfile_4.open ("resultant_lagrangian_strain.txt");
            
    //         for(unsigned int n=0; n<time.get_timestep(); n++)
    //         {
    //             myfile_1<< resultant_cauchy_stress[n]<<std::endl;
    //             myfile_2<< resultant_first_piola_stress[n]<<std::endl;
    //             myfile_3<< resultant_second_piola_stress[n]<<std::endl;
    //             myfile_4<< resultant_lagrangian_strain[n]<<std::endl;
    //         }
    //         myfile_1.close();
    //         myfile_2.close();
    //         myfile_3.close();
    //         myfile_4.close();
    //     }
    // }

    template <int dim>
    unsigned int Solid<dim>::solve()
    {
        vectorType completely_distributed_solution(locally_owned_dofs, mpi_communicator);

        const double success_tol = /*(system_rhs.l2_norm()>1e-10)? 1.e-10 : */1e-3*system_rhs.l2_norm();
        SolverControl solver_control(dof_handler.n_dofs(), success_tol);

        // TrilinosWrappers::SolverBicgstab solver(solver_control);
        TrilinosWrappers::SolverCG solver(solver_control);

        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize(tangent_matrix);
        solver.solve(tangent_matrix, completely_distributed_solution, system_rhs, preconditioner);

        constraints.distribute(completely_distributed_solution);

        solution_update = completely_distributed_solution;

        return solver_control.last_step();
    }

    template <int dim>
    void Solid<dim>::solve_c()
    {
        assemble_system_c();
        vectorType   solution_update_c1(locally_owned_dofs_c, mpi_communicator);
        vectorType   solution_update_c2(locally_owned_dofs_c, mpi_communicator);
        vectorType   solution_update_c3(locally_owned_dofs_c, mpi_communicator);

        vectorType   temp_solution_c1(locally_owned_dofs_c, mpi_communicator);
        vectorType   temp_solution_c2(locally_owned_dofs_c, mpi_communicator);
        vectorType   temp_solution_c3(locally_owned_dofs_c, mpi_communicator);

        temp_solution_c1 = solution_c1;
        temp_solution_c2 = solution_c2;
        temp_solution_c3 = solution_c3;


        const double success_tol_c1 = (system_rhs_c1.l2_norm()==0)? 1.e-3 : 1e-3*system_rhs_c1.l2_norm();
        const double success_tol_c2 = (system_rhs_c2.l2_norm()==0)? 1.e-3 : 1e-3*system_rhs_c2.l2_norm();
        const double success_tol_c3 = (system_rhs_c3.l2_norm()==0)? 1.e-3 : 1e-3*system_rhs_c3.l2_norm();

        SolverControl solver_control_c1 (dof_handler_c.n_dofs(), success_tol_c1);
        SolverControl solver_control_c2 (dof_handler_c.n_dofs(), success_tol_c2);
        SolverControl solver_control_c3 (dof_handler_c.n_dofs(), success_tol_c3);

        // TrilinosWrappers::SolverBicgstab solver_c1 (solver_control_c1);
        // TrilinosWrappers::SolverBicgstab solver_c2 (solver_control_c2);
        // TrilinosWrappers::SolverBicgstab solver_c3 (solver_control_c3);
        TrilinosWrappers::SolverCG solver_c1 (solver_control_c1);
        TrilinosWrappers::SolverCG solver_c2 (solver_control_c2);
        TrilinosWrappers::SolverCG solver_c3 (solver_control_c3);

        TrilinosWrappers::PreconditionAMG preconditioner;
        preconditioner.initialize(mass_matrix);
 
        // std::ostream sout;
        // mass_matrix.print(std::cout, false);

        solver_c1.solve (mass_matrix, solution_update_c1, system_rhs_c1, preconditioner);
        solver_c2.solve (mass_matrix, solution_update_c2, system_rhs_c2, preconditioner);
        solver_c3.solve (mass_matrix, solution_update_c3, system_rhs_c3, preconditioner);

        constraints_c.distribute (solution_update_c1);
        constraints_c.distribute (solution_update_c2);
        constraints_c.distribute (solution_update_c3);

        // temp_solution_c1 += solution_update_c1;
        // temp_solution_c2 += solution_update_c2;
        // temp_solution_c3 += solution_update_c3;

        // double avg_c0 = 0.0;


        IndexSet::ElementIterator it = locally_owned_dofs_c.begin();
        unsigned int n_dofs_per_core = locally_owned_dofs_c.n_elements();

        for(unsigned int i = 0; i < n_dofs_per_core; i++)
        {
            std::vector<double> dc_vec  = {solution_update_c1(*it), solution_update_c2(*it), solution_update_c3(*it)};
            std::vector<double> c_vec   = {temp_solution_c1(*it), temp_solution_c2(*it), temp_solution_c3(*it)};
            std::vector<double> dc_vec_old  = {solution_update_c1(*it), solution_update_c2(*it), solution_update_c3(*it)};


            // if ((1 - (c_vec[0]+c_vec[1]+c_vec[2]) - (dc_vec[0]+dc_vec[1]+dc_vec[2])) < 0)
            // {
            //     dc_vec[0] = dc_vec[0] * (1 - (c_vec[0]+c_vec[1]+c_vec[2]))/abs(dc_vec_old[0]+dc_vec[1]+dc_vec[2]);
            //     dc_vec[1] = dc_vec[1] * (1 - (c_vec[0]+c_vec[1]+c_vec[2]))/abs(dc_vec_old[0]+dc_vec[1]+dc_vec[2]);
            //     dc_vec[2] = dc_vec[2] * (1 - (c_vec[0]+c_vec[1]+c_vec[2]))/abs(dc_vec_old[0]+dc_vec[1]+dc_vec[2]);
            // }

            // AssertThrow((dc_vec[0]+dc_vec[1]+dc_vec[2]) >= 0, "Reverse PT happening");

            // for (unsigned int p = 0; p < 3; p++)
            // {
            //     dc_vec[i]       = trunc(dc_vec[i]);
            //     c_vec[i]        = trunc(c_vec[i]);
            //     dc_vec_old[i]   = trunc(dc_vec_old[i]);
            // }

            // if ((dc_vec[0]+dc_vec[1]+dc_vec[2])>=0)
            {
                if ((c_vec[0]+dc_vec[0]>=0)&&(c_vec[1]+dc_vec[1]>=0)&&(c_vec[2]+dc_vec[2]>=0))
                {
                    if ((c_vec[0]+dc_vec[0]+c_vec[1]+dc_vec[1]+c_vec[2]+dc_vec[2])<=1)
                    {
                        dc_vec[0] = dc_vec[0];
                        dc_vec[1] = dc_vec[1];
                        dc_vec[2] = dc_vec[2];
                    }
                    else
                    {
                        dc_vec[0] = dc_vec[0]*abs((1-(c_vec[0]+c_vec[1]+c_vec[2]))/(abs(dc_vec_old[0])+abs(dc_vec_old[1])+abs(dc_vec_old[2])));
                        dc_vec[1] = dc_vec[1]*abs((1-(c_vec[0]+c_vec[1]+c_vec[2]))/(abs(dc_vec_old[0])+abs(dc_vec_old[1])+abs(dc_vec_old[2])));
                        dc_vec[2] = dc_vec[2]*abs((1-(c_vec[0]+c_vec[1]+c_vec[2]))/(abs(dc_vec_old[0])+abs(dc_vec_old[1])+abs(dc_vec_old[2])));
                    }
                }
                else if ((c_vec[0]+dc_vec[0]<0)&&(c_vec[1]+dc_vec[1]>=0)&&(c_vec[2]+dc_vec[2]>=0))
                {
                    dc_vec[0] = -c_vec[0];
                    if ((c_vec[1]+dc_vec[1]+c_vec[2]+dc_vec[2])<=1)
                    {
                        dc_vec[1] = dc_vec[1];
                        dc_vec[2] = dc_vec[2];
                    }
                    else
                    {
                        dc_vec[1] = dc_vec[1]*abs((1-(c_vec[1]+c_vec[2]))/(abs(dc_vec_old[1])+abs(dc_vec_old[2])));
                        dc_vec[2] = dc_vec[2]*abs((1-(c_vec[1]+c_vec[2]))/(abs(dc_vec_old[1])+abs(dc_vec_old[2])));
                    }
                }
                else if ((c_vec[0]+dc_vec[0]>=0)&&(c_vec[1]+dc_vec[1]<0)&&(c_vec[2]+dc_vec[2]>=0))
                {
                    dc_vec[1] = -c_vec[1];
                    if ((c_vec[0]+dc_vec[0]+c_vec[2]+dc_vec[2])<=1)
                    {
                        dc_vec[0] = dc_vec[0];
                        dc_vec[2] = dc_vec[2];
                    }
                    else
                    {
                        dc_vec[0] = dc_vec[0]*abs((1-(c_vec[0]+c_vec[2]))/(abs(dc_vec_old[0])+abs(dc_vec_old[2])));
                        dc_vec[2] = dc_vec[2]*abs((1-(c_vec[0]+c_vec[2]))/(abs(dc_vec_old[0])+abs(dc_vec_old[2])));
                    }
                }
                else if ((c_vec[0]+dc_vec[0]>=0)&&(c_vec[1]+dc_vec[1]>=0)&&(c_vec[2]+dc_vec[2]<0))
                {
                    dc_vec[2] = -c_vec[2];
                    if ((c_vec[0]+dc_vec[0]+c_vec[1]+dc_vec[1])<=1)
                    {
                        dc_vec[0] = dc_vec[0];
                        dc_vec[1] = dc_vec[1];
                    }
                    else
                    {
                        dc_vec[0] = dc_vec[0]*abs((1-(c_vec[0]+c_vec[1]))/(abs(dc_vec_old[0])+abs(dc_vec_old[1])));
                        dc_vec[1] = dc_vec[1]*abs((1-(c_vec[0]+c_vec[1]))/(abs(dc_vec_old[0])+abs(dc_vec_old[1])));
                    }
                }
                else if ((c_vec[0]+dc_vec[0]<0)&&(c_vec[1]+dc_vec[1]<0)&&(c_vec[2]+dc_vec[2]>=0))
                {
                    dc_vec[0] = -c_vec[0];
                    dc_vec[1] = -c_vec[1];
                    if ((c_vec[2]+dc_vec[2])<=1)
                    {
                        dc_vec[2] = dc_vec[2];
                    }
                    else
                    {
                        dc_vec[2] = 1-c_vec[2];
                    }
                }
                else if ((c_vec[0]+dc_vec[0]<0)&&(c_vec[1]+dc_vec[1]>=0)&&(c_vec[2]+dc_vec[2]<0))
                {
                    dc_vec[0] = -c_vec[0];
                    dc_vec[2] = -c_vec[2];
                    if ((c_vec[1]+dc_vec[1])<=1)
                    {
                        dc_vec[1] = dc_vec[1];
                    }
                    else
                    {
                        dc_vec[1] = 1-c_vec[1];
                    }
                }
                else if ((c_vec[0]+dc_vec[0]>=0)&&(c_vec[1]+dc_vec[1]<0)&&(c_vec[2]+dc_vec[2]<0))
                {
                    dc_vec[1] = -c_vec[1];
                    dc_vec[2] = -c_vec[2];
                    if ((c_vec[0]+dc_vec[0]<=1))
                    {
                        dc_vec[0] = dc_vec[0];
                    }
                    else
                    {
                        dc_vec[0] = 1-c_vec[0];
                    }
                }
                
                else if ((c_vec[0]+dc_vec[0]<0)&&(c_vec[1]+dc_vec[1]<0)&&(c_vec[2]+dc_vec[2]<0))
                {
                    dc_vec[0] = -c_vec[0];
                    dc_vec[1] = -c_vec[1];
                    dc_vec[2] = -c_vec[2];
                }
            }

            trunc(dc_vec[0]);
            trunc(dc_vec[1]);
            trunc(dc_vec[2]);
            trunc(c_vec[0]);
            trunc(c_vec[1]);
            trunc(c_vec[2]);

            // if (dc_vec[0]+dc_vec[1]+dc_vec[2]>=0)
            // {
            temp_solution_c1(*it) = c_vec[0] + dc_vec[0];
            temp_solution_c2(*it) = c_vec[1] + dc_vec[1];
            temp_solution_c3(*it) = c_vec[2] + dc_vec[2];
            // }
            // else
            // {
            //     temp_solution_c1(*it) = c_vec[0];
            //     temp_solution_c2(*it) = c_vec[1];
            //     temp_solution_c3(*it) = c_vec[2];
            // }

            double excess = temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it) - 1;
            if (temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it) > 1)
            {
                if (excess > 0)
                {
                    temp_solution_c1(*it) = temp_solution_c1(*it) - excess/3;
                    temp_solution_c2(*it) = c_vec[1] - excess/3;
                    temp_solution_c3(*it) = c_vec[2] - excess/3;
                }
            }

            temp_solution_c1(*it) = (abs(temp_solution_c1(*it)) < 1e-6)? 0 : trunc(temp_solution_c1(*it));
            temp_solution_c2(*it) = (abs(temp_solution_c2(*it)) < 1e-6)? 0 : trunc(temp_solution_c2(*it));
            temp_solution_c3(*it) = (abs(temp_solution_c3(*it)) < 1e-6)? 0 : trunc(temp_solution_c3(*it));

            // pcout<<"c1:"<<temp_solution_c1(*it)<<std::endl;
            // pcout<<"c2:"<<temp_solution_c2(*it)<<std::endl;
            // pcout<<"c3:"<<temp_solution_c3(*it)<<std::endl;
            // pcout<<"Total:"<<(temp_solution_c1(*it)+temp_solution_c2(*it)+temp_solution_c3(*it))<<std::endl;

            // if (temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it) > 1)
            // {
            //     int p = 0;
            // }

            AssertThrow(temp_solution_c1(*it) >= 0, ExcMessage("c1 < 0"));
            AssertThrow(temp_solution_c2(*it) >= 0, ExcMessage("c2 < 0"));
            AssertThrow(temp_solution_c3(*it) >= 0, ExcMessage("c3 < 0"));
            AssertThrow(temp_solution_c1(*it) <= 1, ExcMessage("c1 > 1"));
            AssertThrow(temp_solution_c2(*it) <= 1, ExcMessage("c2 > 1"));
            AssertThrow(temp_solution_c3(*it) <= 1, ExcMessage("c3 > 1"));

            AssertThrow(temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it) >= 0, ExcMessage("c0 > 1"));
            AssertThrow((temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it)) -1 <= 1e-5, ExcMessage("c0 < 0"));

            // avg_c0 += (1 - (temp_solution_c1(*it) + temp_solution_c2(*it) + temp_solution_c3(*it)))/27;

            it++;
        }

        solution_c1= temp_solution_c1;
        solution_c2= temp_solution_c2;
        solution_c3= temp_solution_c3;

        update_qph_incremental();

        vectorType temp_solution_c0(locally_owned_dofs_c, mpi_communicator);
        temp_solution_c0 = 1;
        temp_solution_c0 -= (temp_solution_c1 + temp_solution_c2 + temp_solution_c3);
        temp_solution_c0.compress(VectorOperation::add);
        solution_c0 = temp_solution_c0;

        // TrilinosScalar total_c0;
        // vectorType zero, one;
        // zero.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        // one.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
        // IndexSet::ElementIterator iter = locally_owned_dofs_c.begin();
        // for(unsigned int i = 0; i < n_dofs_per_core; i++)
        // {
        //     one(*iter) = 1;
        //     iter++;
        // }
        
        // solution_c0.add_and_dot(total_c0, zero, one);
        // pcout   << avg_c0 << std::endl;
        // std::vector<double> average_c0;


        // if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        // {
        //     // average_c0[time.get_timestep()] = avg_c0;


        //     myfile_5.open ("average_c0.txt");
            
        //     for(unsigned int n=0; n<time.get_timestep(); n++)
        //     {
        //         myfile_5<< avg_c0<<std::endl;
        //     }
        //     myfile_5.close();
        // }

         pcout  << "Solving for Volume fractions: \n"
                << solver_control_c1.last_step()<<"+"
                << solver_control_c2.last_step()<<"+"
                << solver_control_c3.last_step()
                << "   CG Solver iterations for C1, C2 and C3."
                << std::endl;
    }

    template <int dim>
    void Solid<dim>::solve_nonlinear_timestep()
    {
        double initial_rhs_norm = 0.;
        unsigned int newton_iteration = 0;
        unsigned int n_iterations=0;
        vectorType  temp_solution_update(locally_owned_dofs, mpi_communicator);
        vectorType  tmp(locally_owned_dofs, mpi_communicator);
        tmp = solution;
        for (; newton_iteration < 3000;   ++newton_iteration)
        {
            make_constraints(newton_iteration);
            assemble_system();

            if (newton_iteration == 0)
            {
                initial_rhs_norm = system_rhs.l2_norm();
                pcout << " Solving for Displacement:   " << std::endl;
            }
            pcout<<"   rhs_norm : "<<system_rhs.l2_norm();

            n_iterations = solve ();
            pcout << "    Number of CG iterations: " << n_iterations<< std::endl;

            temp_solution_update = solution_update;

            tmp += temp_solution_update;
            solution = tmp;

            update_qph_incremental();

            if (newton_iteration > 3 && system_rhs.l2_norm() >= 1)
            {
                pcout << "Diverging... " << std::endl;
                // rhs_norm.clear();

                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            if (newton_iteration > 0 && ((system_rhs.l2_norm() <= 1e-2 * initial_rhs_norm)||(system_rhs.l2_norm() <= 1e-09)))
            {
                pcout << "CONVERGED! " << std::endl;
                break;
            }
            AssertThrow(newton_iteration < 2999, ExcMessage("No convergence in nonlinear solver!"));
        }
    }

    template <int dim>
    void Solid<dim>::run()
    {
        // if (parameters.restart)
        // {
        //     load_checkpoint();
        // }
        // else
        // {
        make_grid(); // generates the geometry and mesh
        // }

        setup_system(); // sets up the system matrices and RHS

        // Applying initial condition for volume fraction c
        vectorType  tmp_solution_c1(locally_owned_dofs_c, mpi_communicator);
        vectorType  tmp_solution_c2(locally_owned_dofs_c, mpi_communicator);
        vectorType  tmp_solution_c3(locally_owned_dofs_c, mpi_communicator);

        if (!parameters.restart)
        {
            // VectorTools::interpolate(dof_handler_c, InitialValues<dim>(1,0), tmp_solution_c1); //initial c
            // VectorTools::interpolate(dof_handler_c, InitialValues<dim>(2,0), tmp_solution_c2); //initial c
            // VectorTools::interpolate(dof_handler_c, InitialValues<dim>(3,0), tmp_solution_c3); //initial c
            solution_c1= tmp_solution_c1;
            solution_c2= tmp_solution_c2;
            solution_c3= tmp_solution_c3;

            update_qph_incremental();

            solve_nonlinear_timestep();
            // output_resultant_stress();
            output_results();
        }
        else
        {
            load_time();
            
            for (int t = 0; t<checkpoint_timestep; ++t)
            {
                time.increment();
            }
        }

        time.increment();

        // computed actual time integration to update displacement and c
        while (time.current() <= time.end() )
        {
            pcout   << std::endl
                    << "Time step #" << time.get_timestep() << "; "
                    << "advancing to t = " << time.current() << "."
                    << std::endl;

            solve_c();

            solve_nonlinear_timestep();

            if((time.get_timestep()%25 == 0 && time.get_timestep()<=100)||(time.get_timestep()%200 == 0 && time.get_timestep()>=100))
            {
                output_results();
            }
            // output_resultant_stress();
            if(time.get_timestep()%1000==0)
            {
                save_checkpoint();
            }

            time.increment();
            // int result = time.get_timestep();
        }
    }
}//namespace PhaseField

// // Defining a global instance of input parameters
// PhaseField::Parameters::AllParameters parameters;

// // Defining a global variable for rotation
// const dealii::Tensor<2, 3, double> Rot_mat()
// {
//     static dealii::Tensor<2, 3, double> Rot_mat = dealii::Physics::Transformations::Rotations::rotation_matrix_3d(dealii::Point<3, double>(0.0, 0.0, 1.0), -1.0472);
//     return Rot_mat;
// }



/////////////////////
int main (int argc, char *argv[])
{
    try
    {
        using namespace dealii;
        using namespace PhaseField;

        Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

        std::string param_file;
        param_file = "parameters.prm";

        parameters.read_prm(param_file);

        Solid<3> solid_3d;
        solid_3d.run();
    }
    catch (std::exception &exc)
    {
        std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
        std::cerr << "Exception on processing: " << std::endl << exc.what()
                    << std::endl << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;

        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
        std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                    << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
        return 1;
    }

    return 0;
}
