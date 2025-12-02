/**
 * @file solid.h
 * @brief Main solver class for coupled phase-field and mechanical simulation
 *
 * This header defines the Solid class which manages the entire simulation
 * workflow:
 * - Mesh generation and refinement
 * - System assembly for both mechanical and phase-field equations
 * - Nonlinear Newton-Raphson solver for mechanical equilibrium
 * - Phase-field evolution solver
 * - Quadrature point history management
 * - Output generation
 *
 * The class uses parallel distributed computing with MPI and Trilinos sparse
 * linear algebra. Template implementation is in ../src/solid.C
 */
#ifndef SOLID_CORE_H
#define SOLID_CORE_H

#include <fstream>
#include <iostream>
// #include <string>

#include "boundarydisplacement.h"
#include "dealiiheaders.h"
#include "fesystem_str.h"
#include "pointhistory.h"
#include "timestep.h"

namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim> class Solid {
public:
  /**
   * @brief Constructor
   * @param input_file Path to parameter file (typically parameters.prm)
   */
  Solid(const std::string &input_file);

  /**
   * @brief Destructor
   */
  virtual ~Solid();

  /**
   * @brief Main simulation driver
   *
   * Orchestrates the entire simulation workflow:
   * - Grid generation
   * - System setup
   * - Time stepping loop
   * - Solving coupled equations
   * - Output generation
   */
  void run();

private:
  /** @brief Generate and refine the computational mesh */
  void make_grid();

  /** @brief Setup DOF handlers, constraints, and system matrices */
  void system_setup();

  /** @brief Assemble the mechanical system matrix and RHS */
  void assemble_system();

  /** @brief Apply boundary conditions and constraints
   *  @param it_nr Newton iteration number
   */
  void make_constraints(const int &it_nr);

  /** @brief Solve nonlinear mechanical problem using Newton-Raphson */
  void solve_nonlinear_timestep();

  /** @brief Solve linear system for displacement increment
   *  @return Number of solver iterations
   */
  unsigned int solve();

  /** @brief Assemble system for phase-field evolution equations */
  void assemble_system_c();

  /** @brief Solve phase-field evolution equations */
  void solve_c();

  /** @brief Setup quadrature point history data structure */
  void setup_qph();

  /** @brief Update material state at all quadrature points */
  void update_qph_incremental();

  /** @brief Write VTU output files for visualization */
  void output_results() const;

  /** @brief Compute and output resultant stresses on boundaries */
  void output_resultant_stress();

  /** @brief Write quadrature point data to files */
  void output_quad();

  /** @brief Write quadrature data for specific increment
   *  @param _currentIncrement Current time step number
   */
  void writeQuadratureOutput(unsigned int _currentIncrement);

  // Parameters::AllParameters        parameters;
  //
  // Time                             time;  // variable of type class 'Time'
  // //TimerOutput                      timer;

  MPI_Comm mpi_communicator;
  parallel::distributed::Triangulation<dim> triangulation;
  // ConditionalOStream               pcout;

  Parameters::AllParameters parameters;

  Time time; // variable of type class 'Time'
  ConditionalOStream pcout;
  mutable TimerOutput timer;

  const unsigned int degree;   // degree of polynomial of shape functions
  const FESystem<dim> fe;      // fe object
  DoFHandler<dim> dof_handler; // we have used two dof_handler: one for
                               // mechanics another for order parameter
  const unsigned int
      dofs_per_cell; // no of dofs per cell for the mechanics problem
  const FEValuesExtractors::Vector u_fe;
  const QGauss<dim> qf_cell;             // quadrature points in the cell
  const QGauss<dim - 1> qf_face;         // quadrature points at the face
  const unsigned int n_q_points;         // no of quadrature points in the cell
  const unsigned int n_q_points_f;       // no of quadrature points at the face
  AffineConstraints<double> constraints; // constraint object

  // FE_DGQ<dim>                     history_fe;
  // DoFHandler<dim>                    history_dof_handler;

  std::vector<PointHistory<dim>> quadrature_point_history;

  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  matrixType tangent_matrix; // tangent stiffenss matrix
  vectorType
      system_rhs; // system right hand side or residual of mechanics problem
  vectorType solution;        // solution vector for displacement
  vectorType solution_update; // another vector containing the displacement soln

  const unsigned int degree_c;        // degree of polynomial for c
  FE_Q<dim> fe_c;                     // fe object for c
  DoFHandler<dim> dof_handler_c;      // another dof_handler for c
  const unsigned int dofs_per_cell_c; // dof per c cell
  const QGauss<dim> qf_cell_c;
  const unsigned int n_q_points_c;
  AffineConstraints<double> constraints_c;
  IndexSet locally_owned_dofs_c;
  IndexSet locally_relevant_dofs_c;
  // const std::vector<DynamicSparsityPattern::size_type>

  DoFHandler<dim> history_dof_handler;
  FE_DGQ<dim> history_fe;

  matrixType mass_matrix;
  vectorType system_rhs_c1, system_rhs_c2, system_rhs_c3;
  vectorType solution_c0, solution_c1, solution_c2, solution_c3;
  vectorType old_solution_c1, old_solution_c2, old_solution_c3;
  vectorType solution_update_c1, solution_update_c2, solution_update_c3;

  Vector<double> resultant_cauchy_stress;
  Vector<double> resultant_first_piola_stress;
  Vector<double> resultant_second_piola_stress;
  Vector<double> static_cauchy_stress;
  Vector<double> static_first_piola_stress;
  Vector<double> static_second_piola_stress;
  Vector<double> resultant_lagrangian_strain;
  Vector<double> static_lagrangian_strain;
  Vector<double> order_parameter;
  Vector<double> static_order_parameter;
  bool apply_strain;
  double load_step;
  double load;

  /** @brief Directory for output files */
  std::string output_directory;
  bool suppress_file_output;

  std::vector<std::vector<double>> outputQuadrature;
};

} // namespace PhaseField
#include "../src/solid.C"
#endif
