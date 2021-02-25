#ifndef SOLID_CORE_H
#define SOLID_CORE_H

#include <iostream>
#include <fstream>

#include "dealiiheaders.h"
#include "fesystem_str.h"
#include "timestep.h"
#include "pointhistory.h"
#include "boundarydisplacement.h"


namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim>
class Solid
{
public:
  Solid(const std::string &input_file);

  virtual
  ~Solid();

  void
  run();

private:

  void    make_grid();
  void    system_setup();
  void    assemble_system();
  void    make_constraints(const int &it_nr);
  void    solve_nonlinear_timestep();
  unsigned int    solve();
  void    assemble_system_c();
  void    solve_c();
  void    setup_qph();
  void    update_qph_incremental();
  void    output_results() const;
  void    output_resultant_stress();


  // Parameters::AllParameters        parameters;
  //
  // Time                             time;  // variable of type class 'Time'
  // //TimerOutput                      timer;


  MPI_Comm                         mpi_communicator;
  parallel::distributed::Triangulation<dim> triangulation;
  // ConditionalOStream               pcout;

  Parameters::AllParameters        parameters;

  Time                             time;  // variable of type class 'Time'
  //TimerOutput                      timer;

  ConditionalOStream               pcout;

  const unsigned int               degree; // degree of polynomial of shape functions
  const FESystem<dim>              fe; // fe object
  DoFHandler<dim>                 dof_handler; // we have used two dof_handler: one for mechanics another for order parameter
  const unsigned int               dofs_per_cell;   // no of dofs per cell for the mechanics problem
  const FEValuesExtractors::Vector   u_fe;
  const QGauss<dim>                qf_cell;  // quadrature points in the cell
  const QGauss<dim - 1>            qf_face;  // quadrature points at the face
  const unsigned int               n_q_points;  // no of quadrature points in the cell
  const unsigned int               n_q_points_f; // no of quadrature points at the face
  AffineConstraints<double>                constraints;  // constraint object

  // FE_DGQ<dim>                     history_fe;
  // DoFHandler<dim>                    history_dof_handler;

  std::vector<PointHistory<dim> >  quadrature_point_history;

  IndexSet                         locally_owned_dofs;
  IndexSet                         locally_relevant_dofs;


  matrixType                   tangent_matrix;  // tangent stiffenss matrix
  vectorType                system_rhs;  // system right hand side or residual of mechanics problem
  vectorType                solution;  // solution vector for displacement
  vectorType                solution_update; // another vector containing the displacement soln


  const unsigned int               degree_c; // degree of polynomial for c
  FE_Q<dim>                        fe_c;  // fe object for c
  DoFHandler<dim>                  dof_handler_c; //another dof_handler for c
  const unsigned int               dofs_per_cell_c; // dof per c cell
  const QGauss<dim>                qf_cell_c;
  const unsigned int               n_q_points_c;
  AffineConstraints<double>        constraints_c;
  IndexSet                         locally_owned_dofs_c;
  IndexSet                         locally_relevant_dofs_c;
  // const std::vector<DynamicSparsityPattern::size_type>

  DoFHandler<dim>                    history_dof_handler;
  FE_DGQ<dim>                     history_fe;

  matrixType                  mass_matrix;
  vectorType                  system_rhs_c1, system_rhs_c2, system_rhs_c3;
  vectorType                  solution_c0, solution_c1, solution_c2, solution_c3;
  vectorType                  old_solution_c1, old_solution_c2, old_solution_c3;
  vectorType                  solution_update_c1, solution_update_c2, solution_update_c3 ;


  Vector<double>                   resultant_cauchy_stress;
  Vector<double>                   resultant_first_piola_stress;
  Vector<double>                   resultant_second_piola_stress;
  Vector<double>                   static_cauchy_stress;
  Vector<double>                   static_first_piola_stress;
  Vector<double>                   static_second_piola_stress;
  Vector<double>                   resultant_lagrangian_strain;
  Vector<double>                   static_lagrangian_strain;
  Vector<double>                   order_parameter;
  Vector<double>                   static_order_parameter;
  bool                             apply_strain;
  double                           load_step;
  double                           load;

};

template <int dim>
Solid<dim>::Solid(const std::string &input_file)
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator,
                 typename Triangulation<dim>::MeshSmoothing
                 (Triangulation<dim>::smoothing_on_refinement |
                  Triangulation<dim>::smoothing_on_coarsening)),

  parameters(input_file),
  time(parameters.end_time, parameters.delta_t),

  pcout (std::cout,
        (Utilities::MPI::this_mpi_process(mpi_communicator)
         == 0)),
//    timer(mpi_communicator,
//          pcout,
//          TimerOutput::summary,
//          TimerOutput::wall_times),

  degree(parameters.poly_degree),
  fe(FE_Q<dim>(parameters.poly_degree), dim), // displacement
  dof_handler(triangulation),
  dofs_per_cell (fe.dofs_per_cell),
  u_fe(0),
  qf_cell(parameters.quad_order),
  qf_face(parameters.quad_order),
  n_q_points (qf_cell.size()),
  n_q_points_f (qf_face.size()),

  degree_c(parameters.poly_degree),
  fe_c (parameters.poly_degree),
  dof_handler_c (triangulation),
  dofs_per_cell_c(fe_c.dofs_per_cell),
  qf_cell_c(parameters.quad_order),
  n_q_points_c (qf_cell_c.size()),

  history_dof_handler (triangulation),
  history_fe (parameters.poly_degree),
  apply_strain(false),
  load_step(1),
  load(0.0)
{}

//destructor
template <int dim>
Solid<dim>::~Solid()
{
  dof_handler.clear();
  dof_handler_c.clear();
}

//make_grid
template <int dim>
void Solid<dim>::make_grid()
{

  std::vector< unsigned int > repetitions(dim, 30);
  if (dim == 3)
  repetitions[dim-2] = 10;
  repetitions[dim-3] = 5;

  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            repetitions,
                                            Point<dim>(0.0, 0.0, 0.0),
                                            Point<dim>(0.5, 1.0, 3.0),
                                               true);
}


//system_setup
template <int dim>
void Solid<dim>::system_setup()
{
  dof_handler.distribute_dofs(fe);
  dof_handler_c.distribute_dofs (fe_c);
  history_dof_handler.distribute_dofs (history_fe);

  const unsigned int n_dofs = dof_handler.n_dofs(),
                      n_dofs_c  = dof_handler_c.n_dofs();

       pcout     << "   Number of active cells: "
                     << triangulation.n_active_cells()
                     << std::endl
                     << "   Total number of cells: "
                     << triangulation.n_cells()
                     << std::endl
                        << "   Number of degrees of freedom: "
                     << n_dofs + n_dofs_c
                     << " (" << n_dofs << '+' << n_dofs_c << ')'
                     << std::endl;

  locally_owned_dofs = dof_handler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler  ,  locally_relevant_dofs);

  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);

  {
//      DoFTools::make_periodicity_constraints(dof_handler,
//                                             /*b_id*/ 0,
//                                             /*b_id*/ 1,
//                                             /*direction*/ 0,
//                                             constraints);
//      DoFTools::make_periodicity_constraints(dof_handler,
//                                             /*b_id*/ 2,
//                                             /*b_id*/ 3,
//                                             /*direction*/ 1,
//                                             constraints);
//      DoFTools::make_periodicity_constraints(dof_handler,
//                                               /*b_id*/ 4,
//                                               /*b_id*/ 5,
//                                               /*direction*/ 2,
//                                               constraints);

 }
 constraints.close ();

  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
  Utilities::MPI::all_gather(mpi_communicator, dof_handler.locally_owned_dofs()),

  // SparsityTools::distribute_sparsity_pattern (dsp,
  //                                             dof_handler.n_locally_owned_dofs_per_processor(),
  //                                             mpi_communicator,
  //                                             locally_relevant_dofs);
  SparsityTools::distribute_sparsity_pattern (dsp,
                                              dof_handler.locally_owned_dofs(),
                                              mpi_communicator,
                                              locally_relevant_dofs);
  // Utilities::MPI::all_gather(mpi_communicator, locally_owned_dofs);

  tangent_matrix.reinit(locally_owned_dofs,
                        locally_owned_dofs,
                        dsp,
                        mpi_communicator);

  solution.reinit(locally_owned_dofs,
                    locally_relevant_dofs,
                    mpi_communicator);
  solution_update.reinit(locally_owned_dofs,
                         locally_relevant_dofs,
                         mpi_communicator);

  system_rhs.reinit(locally_owned_dofs,
                    mpi_communicator);


  locally_owned_dofs_c = dof_handler_c.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dof_handler_c  ,  locally_relevant_dofs_c);

  constraints_c.clear ();
  constraints_c.reinit (locally_relevant_dofs_c);
  {

//        DoFTools::make_periodicity_constraints(dof_handler_c,
//                                               /*b_id*/ 0,
//                                               /*b_id*/ 1,
//                                               /*direction*/ 0,
//                                               constraints_c);
//        DoFTools::make_periodicity_constraints(dof_handler_c,
//                                                 /*b_id*/ 2,
//                                               /*b_id*/ 3,
//                                               /*direction*/ 1,
//                                               constraints_c);
//        DoFTools::make_periodicity_constraints(dof_handler_c,
//                                               /*b_id*/ 4,
//                                               /*b_id*/ 5,
//                                               /*direction*/ 2,
//                                               constraints_c);

  }
  constraints_c.close ();

  DynamicSparsityPattern dsp_c(locally_relevant_dofs_c);
  DoFTools::make_sparsity_pattern (dof_handler_c, dsp_c, constraints_c, true);
  Utilities::MPI::all_gather(mpi_communicator, dof_handler_c.locally_owned_dofs());

  // SparsityTools::distribute_sparsity_pattern (dsp_c,
  //                                             dof_handler_c.n_locally_owned_dofs_per_processor(),
  //                                             mpi_communicator,
  //                                             locally_relevant_dofs_c);
  SparsityTools::distribute_sparsity_pattern (dsp_c,
                                              dof_handler_c.locally_owned_dofs(),
                                              mpi_communicator,
                                              locally_relevant_dofs_c);
  // Utilities::MPI::all_gather(mpi_communicator, locally_owned_dofs_c);

  mass_matrix.reinit (locally_owned_dofs_c, locally_owned_dofs_c, dsp_c, mpi_communicator);


  solution_c0.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
  solution_c1.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
  solution_c2.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);
  solution_c3.reinit(locally_owned_dofs_c, locally_relevant_dofs_c, mpi_communicator);


  system_rhs_c1.reinit(locally_owned_dofs_c, mpi_communicator);
  system_rhs_c2.reinit(locally_owned_dofs_c, mpi_communicator);
  system_rhs_c3.reinit(locally_owned_dofs_c, mpi_communicator);

  resultant_cauchy_stress.reinit(100000);
  static_cauchy_stress.reinit   (100000);
  resultant_first_piola_stress.reinit(100000);
  static_first_piola_stress.reinit   (100000);
  resultant_second_piola_stress.reinit(100000);
  static_second_piola_stress.reinit   (100000);
  resultant_lagrangian_strain.reinit(100000);
  static_lagrangian_strain.reinit(100000);
  order_parameter.reinit(100000);
  static_order_parameter.reinit(100000);

 setup_qph();


}


//make_constraints
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

// // Fixing points or lines

//        const double tol_boundary = 0.01;
//        typename DoFHandler<dim>::active_cell_iterator
//         cell = dof_handler.begin_active(),
//        endc = dof_handler.end();
//        for (; cell!=endc; ++cell)
//            if (cell->is_locally_owned())
//              {
//                  for (unsigned int v=0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
//
//                  if     ((std::abs(cell->vertex(v)[0] - 0.5) < tol_boundary) &&
//                          (std::abs(cell->vertex(v)[1] - 0.5) < tol_boundary) &&
//                          (std::abs(cell->vertex(v)[2] - 1.0) < tol_boundary))
//                       {
//                            constraints.add_line(cell->vertex_dof_index(v, 0));
//                            constraints.add_line(cell->vertex_dof_index(v, 1));
//                            constraints.add_line(cell->vertex_dof_index(v, 2));
//                       }
//            else if     (/*(std::abs(cell->vertex(v)[0] - 0.0) < tol_boundary) &&*/
//                         (std::abs(cell->vertex(v)[1] - 0.0) < tol_boundary) &&
//                         (std::abs(cell->vertex(v)[2] - 0.0) < tol_boundary))
//                      {
//                         // constraints.add_line(cell->vertex_dof_index(v, 0));
//                         constraints.add_line(cell->vertex_dof_index(v, 1));
//                         //constraints.add_line(cell->vertex_dof_index(v, 2));
//                      }
//
//           else if      ((std::abs(cell->vertex(v)[0] - 0.0) < tol_boundary) &&
//                         /*(std::abs(cell->vertex(v)[1] - 0.0) < tol_boundary) &&*/
//                         (std::abs(cell->vertex(v)[2] - 0.0) < tol_boundary))
//                      {
//                         constraints.add_line(cell->vertex_dof_index(v, 0));
//                        //constraints.add_line(cell->vertex_dof_index(v, 1));
//                        //constraints.add_line(cell->vertex_dof_index(v, 2));
//                     }
//           else if      ((std::abs(cell->vertex(v)[0] - 0.0) < tol_boundary) &&
//                         (std::abs(cell->vertex(v)[1] - 1.0) < tol_boundary) &&
//                         (std::abs(cell->vertex(v)[2] - 0.0) < tol_boundary))
//                     {
//                        constraints.add_line(cell->vertex_dof_index(v, 0));
//                         //constraints.add_line(cell->vertex_dof_index(v, 1));
//                        //constraints.add_line(cell->vertex_dof_index(v, 2));
//                     }
//              }

// Fixing external surfaces
                  {
                       const int boundary_id = 0;

                  if (apply_dirichlet_bc == true)
                    VectorTools::interpolate_boundary_values(dof_handler,
                                                             boundary_id,
                                                             ZeroFunction<dim>(dim),
                                                             constraints,
                                                             (fe.component_mask(x_displacement)));
                  else
                    VectorTools::interpolate_boundary_values(dof_handler,
                                                             boundary_id,
                                                             ZeroFunction<dim>(dim),
                                                             constraints,
                                                             (fe.component_mask(x_displacement)));
                  }
//
//                    {
//                        const int boundary_id = 1;
//
//                         if (apply_dirichlet_bc == true)
//                          VectorTools::interpolate_boundary_values(dof_handler,
//                                                                   boundary_id,
//                                                                      BoundaryDisplacement<dim>(0, timestep),
//                                                                   constraints,
//                                                                   fe.component_mask(x_displacement));
//                        else
//                          VectorTools::interpolate_boundary_values(dof_handler,
//                                                                   boundary_id,
//                                                                   ZeroFunction<dim>(dim),
//                                                                   constraints,
//                                                                   fe.component_mask(x_displacement));
//                   }

//                    {
//                          const int boundary_id = 2;
//
//                     if (apply_dirichlet_bc == true)
//                       VectorTools::interpolate_boundary_values(dof_handler,
//                                                                boundary_id,
//                                                                ZeroFunction<dim>(dim),
//                                                               constraints,
//                                                                fe.component_mask(y_displacement));
//                     else
//                       VectorTools::interpolate_boundary_values(dof_handler,
//                                                                boundary_id,
//                                                                ZeroFunction<dim>(dim),
//                                                                constraints,
//                                                                fe.component_mask(y_displacement));
//                     }

//                    {
//                           const int boundary_id = 3;
//
//                      if (apply_dirichlet_bc == true)
//                        VectorTools::interpolate_boundary_values(dof_handler,
//                                                                 boundary_id,
//                                                               BoundaryDisplacement<dim>(1, timestep),
//                                                                constraints,
//                                                                 fe.component_mask(y_displacement));
//                      else
//                        VectorTools::interpolate_boundary_values(dof_handler,
//                                                                 boundary_id,
//                                                                 ZeroFunction<dim>(dim),
//                                                                 constraints,
//                                                                 fe.component_mask(y_displacement));
//                      }

                  {

                      const int boundary_id = 4;

                        if (apply_dirichlet_bc == true)
                          VectorTools::interpolate_boundary_values(dof_handler,
                                                                   boundary_id,
                                                                   ZeroFunction<dim>(dim),
                                                                   constraints,
                                                                   fe.component_mask(x_displacement)
                                                                 |
                                                                 fe.component_mask(y_displacement)
                                                                 |
                                                                 fe.component_mask(z_displacement));
                        else
                          VectorTools::interpolate_boundary_values(dof_handler,
                                                                   boundary_id,
                                                                   ZeroFunction<dim>(dim),
                                                                   constraints,
                                                                 fe.component_mask(x_displacement)
                                                                 |
                                                                 fe.component_mask(y_displacement)
                                                                 |
                                                                 fe.component_mask(z_displacement));
                   }
                  {

                      const int boundary_id = 5;

                        if (apply_dirichlet_bc == true)
                          VectorTools::interpolate_boundary_values(dof_handler,
                                                                   boundary_id,
                                                                 BoundaryDisplacement<dim>(2, timestep),
                                                                   constraints,
                                                                   fe.component_mask(z_displacement));
                        else
                          VectorTools::interpolate_boundary_values(dof_handler,
                                                                   boundary_id,
                                                                   ZeroFunction<dim>(dim),
                                                                   constraints,
                                                                   fe.component_mask(z_displacement));
                   }
                 {

                     const int boundary_id = 5;

                       if (apply_dirichlet_bc == true)
                         VectorTools::interpolate_boundary_values(dof_handler,
                                                                  boundary_id,
                                                                 ZeroFunction<dim>(dim),
                                                                  constraints,
                                                                  fe.component_mask(x_displacement)
                                                                 |
                                                                 fe.component_mask(y_displacement));
                       else
                         VectorTools::interpolate_boundary_values(dof_handler,
                                                                  boundary_id,
                                                                  ZeroFunction<dim>(dim),
                                                                  constraints,
                                                                  fe.component_mask(x_displacement)
                                                                 |
                                                                 fe.component_mask(y_displacement));
                  }

// In the case of periodic conditions these lines shold be uncommented.

//{
//        DoFTools::make_periodicity_constraints(dof_handler,
//                                              /*b_id*/ 0,
//                                              /*b_id*/ 1,
//                                              /*direction*/ 0,
//                                              constraints);
//        DoFTools::make_periodicity_constraints(dof_handler,
//                                              /*b_id*/ 2,
//                                              /*b_id*/ 3,
//                                              /*direction*/ 1,
//                                                 constraints);
//       DoFTools::make_periodicity_constraints(dof_handler,
//                                               /*b_id*/ 4,
//                                               /*b_id*/ 5,
//                                               /*direction*/ 2,
//                                               constraints);
//      }
//     constraints.close();

//     // Applying boundary displacement when having periodic conditions
/*
    {
       IndexSet selected_dofs_x;
       std::set< types::boundary_id > boundary_ids_x= std::set<types::boundary_id>();
               boundary_ids_x.insert(0);

       DoFTools::extract_boundary_dofs(dof_handler,
                                      fe.component_mask(x_displacement),
                                           selected_dofs_x,
                                           boundary_ids_x);
       unsigned int nb_dofs_face_x = selected_dofs_x.n_elements();
       IndexSet::ElementIterator dofs_x = selected_dofs_x.begin();

       double relative_displacement_x;

       if(timestep<80)
       relative_displacement_x = 5e-4;
       else if(timestep%100==0)
       relative_displacement_x = 5e-4;
       else
       relative_displacement_x = 0.0;

       for(unsigned int i = 0; i < nb_dofs_face_x; i++)
       {
        constraints.add_line (*dofs_x);
         constraints.set_inhomogeneity(*dofs_x, (apply_dirichlet_bc ? relative_displacement_x : 0.0));
           dofs_x++;
       }
     }


   {
       IndexSet selected_dofs_y;
       std::set< types::boundary_id > boundary_ids_y= std::set<types::boundary_id>();
               boundary_ids_y.insert(2);

       DoFTools::extract_boundary_dofs(dof_handler,
                                      fe.component_mask(y_displacement),
                                         selected_dofs_y,
                                         boundary_ids_y);
       unsigned int nb_dofs_face_y = selected_dofs_y.n_elements();
       IndexSet::ElementIterator dofs_y = selected_dofs_y.begin();

       double relative_displacement_y = 0.0;

       if(timestep<80)
        relative_displacement_y = 5e-4;
       else if(timestep%100==0)
        relative_displacement_y = 5e-4;
       else
        relative_displacement_y = 0.0;


       for(unsigned int i = 0; i < nb_dofs_face_y; i++)
       {
         constraints.add_line (*dofs_y);
         constraints.set_inhomogeneity(*dofs_y,(apply_dirichlet_bc ? relative_displacement_y : 0.0) );
           dofs_y++;
       }

     }

   {
     IndexSet selected_dofs_z;
     std::set< types::boundary_id > boundary_ids_z= std::set<types::boundary_id>();
             boundary_ids_z.insert(4);

     DoFTools::extract_boundary_dofs(dof_handler,
                                    fe.component_mask(z_displacement),
                                         selected_dofs_z,
                                         boundary_ids_z);
     unsigned int nb_dofs_face_z = selected_dofs_z.n_elements();
     IndexSet::ElementIterator dofs_z = selected_dofs_z.begin();

     double relative_displacement_z;
     relative_displacement_z = 0.0;

     if(timestep<75)
      relative_displacement_z = 10e-4;
     else if(timestep%10==0)
      relative_displacement_z = 10e-4;
     else
      relative_displacement_z = 0.0;


     for(unsigned int i = 0; i < nb_dofs_face_z; i++)
     {
       constraints.add_line (*dofs_z);
         constraints.set_inhomogeneity(*dofs_z,(apply_dirichlet_bc ? relative_displacement_z : 0.0));
         dofs_z++;
     }
  }
*/
  constraints.close();

}


//assemble_system
template <int dim>
void Solid<dim>::assemble_system ()
{
  tangent_matrix = 0;
  system_rhs = 0;

 FEValues<dim> fe_values (fe, qf_cell,
                          update_values   | update_gradients |
                          update_quadrature_points | update_JxW_values);


 FEFaceValues<dim> fe_face_values (fe, qf_face,
                                   update_values         | update_quadrature_points  |
                                   update_normal_vectors | update_JxW_values);

 FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
 Vector<double>       cell_rhs (dofs_per_cell);

 std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

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

     PointHistory<dim> *lqph =
           reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
         {
          const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();
          const Tensor<2, dim> tau   = lqph[q_point].get_tau();
          const SymmetricTensor<2, dim> symm_tau         = lqph[q_point].get_tau();
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

                    cell_matrix(i, j) += symm_grad_Nx[i] * Jc // The material contribution:
                                                * symm_grad_Nx[j] * JxW;
                   if (component_i == component_j) // geometrical stress contribution
                    cell_matrix(i, j) += grad_Nx[i][component_i] * tau
                                                * grad_Nx[j][component_j] * JxW;
            }

                    cell_rhs(i) -= symm_grad_Nx[i] * symm_tau * JxW;
          }
     }




// Applying force on the external faces

//       // if (time.get_timestep()>3)
//       {
//
//           for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
//
//            if (cell->face(face_number)->at_boundary()
//                    &&
//               ( (cell->face(face_number)->boundary_id() == 0 )
//                     ||
//                 (cell->face(face_number)->boundary_id() == 1 )
//                     ||
//                 (cell->face(face_number)->boundary_id() == 2 )
//                     ||
//                 (cell->face(face_number)->boundary_id() == 3 )))
//             {
//               fe_face_values.reinit (cell, face_number);
//
//               for (unsigned int q_point=0; q_point<n_q_points_f; ++q_point)
//                 {
//                   const Tensor<2, dim> F_inv_tr = lqph[q_point].get_F_inv_tr();
//
//                     const double J= lqph[q_point].get_det_F();
//
//                           double normal_stress;
//                         if (load_step<12)
//                              normal_stress= -load_step * J * (F_inv_tr*fe_face_values.normal_vector(q_point)).norm();
//                         else
//                              normal_stress= -12 * J * (F_inv_tr*fe_face_values.normal_vector(q_point)).norm();
//
//                         const Tensor<1, dim> traction  = normal_stress* fe_face_values.normal_vector(q_point);
//
//
//                  for (unsigned int i=0; i<dofs_per_cell; ++i)
//                   {
//                       const unsigned int
//                     component_i = fe.system_to_component_index(i).first;
//                     cell_rhs(i) += (traction[component_i] *
//                                     fe_face_values.shape_value(i,q_point) *
//                                     fe_face_values.JxW(q_point));
//                   }
//                 }
//               }
//
//            if (cell->face(face_number)->at_boundary()
//                     &&
//             (cell->face(face_number)->boundary_id() == 5)
//
//                 )
//                 {
//                   fe_face_values.reinit (cell, face_number);
//                  for (unsigned int q_point=0; q_point<n_q_points_f; ++q_point)
//                     {
//                      const Tensor<2, dim> F_inv_tr = lqph[q_point].get_F_inv_tr();
//
//                         const double J= lqph[q_point].get_det_F();
//
//                         double normal_stress;
//                         if (load_step<4)
//                            normal_stress= -2*load_step * J * (F_inv_tr*fe_face_values.normal_vector(q_point)).norm();
//                         else
//                        normal_stress= -8.0 * J * (F_inv_tr*fe_face_values.normal_vector(q_point)).norm();
//
//                      const Tensor<1, dim> traction  = normal_stress* fe_face_values.normal_vector(q_point);
//
//
//                       for (unsigned int i=0; i<dofs_per_cell; ++i)
//                       {
//                           const unsigned int
//                         component_i = fe.system_to_component_index(i).first;
//                         cell_rhs(i) += (traction[component_i] *
//                                         fe_face_values.shape_value(i,q_point) *
//                                         fe_face_values.JxW(q_point));
//                       }
//                     }
//                   }
//
//
//              }


    cell->get_dof_indices (local_dof_indices);
    constraints.distribute_local_to_global (cell_matrix,
                                            cell_rhs,
                                            local_dof_indices,
                                            tangent_matrix,
                                            system_rhs);

   }

 tangent_matrix.compress (VectorOperation::add);
 system_rhs.compress (VectorOperation::add);

}

// Assemplying system matrix and RHS for the phase-field kinetic equation
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

  FullMatrix<double>   cell_mass_matrix    (dofs_per_cell_c, dofs_per_cell_c);
  Vector<double>   cell_rhs_c1         (dofs_per_cell_c);
  Vector<double>   cell_rhs_c2         (dofs_per_cell_c);
  Vector<double>   cell_rhs_c3         (dofs_per_cell_c);

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell_c);

  std::vector<double> phi (dofs_per_cell_c);

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler_c.begin_active(),
  endc = dof_handler_c.end();
  for (; cell!=endc; ++cell)
    if (cell->is_locally_owned())
    {
      fe_values_c.reinit(cell);

      cell_mass_matrix = 0;
      cell_rhs_c1=0;
      cell_rhs_c2=0;
      cell_rhs_c3=0;

      PointHistory<dim> *lqph =
                    reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

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


      cell->get_dof_indices (local_dof_indices);
      constraints_c.distribute_local_to_global (cell_mass_matrix,
                                                local_dof_indices,
                                                mass_matrix);
      constraints_c.distribute_local_to_global (cell_rhs_c1,
                                                local_dof_indices,
                                                system_rhs_c1);
      constraints_c.distribute_local_to_global (cell_rhs_c2,
                                                local_dof_indices,
                                                system_rhs_c2);
      constraints_c.distribute_local_to_global (cell_rhs_c3,
                                                local_dof_indices,
                                                system_rhs_c3);

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
  for (typename Triangulation<dim>::active_cell_iterator
       cell = triangulation.begin_active();
       cell != triangulation.end(); ++cell)
    if (cell->is_locally_owned())
      ++our_cells;
  triangulation.clear_user_data();
  {
    std::vector<PointHistory<dim> > tmp;
    tmp.swap (quadrature_point_history);
  }
  quadrature_point_history.resize (our_cells * n_q_points);

    unsigned int history_index = 0;
    for (typename Triangulation<dim>::active_cell_iterator
            cell = triangulation.begin_active();
            cell != triangulation.end(); ++cell)
     if (cell->is_locally_owned())
      {
        cell->set_user_pointer(&quadrature_point_history[history_index]);
        history_index += n_q_points;
      }

    Assert(history_index == quadrature_point_history.size(),
           ExcInternalError());
  }

  for (typename Triangulation<dim>::active_cell_iterator
          cell = triangulation.begin_active();
          cell != triangulation.end(); ++cell)
   if (cell->is_locally_owned())
    {
      PointHistory<dim> *lqph =
        reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

      Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
      Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        lqph[q_point].setup_lqp(parameters);
    }
}


//update_qph_incremental
template <int dim>
void Solid<dim>::update_qph_incremental()
{
  FEValues<dim> fe_values (fe, qf_cell,
                               update_values | update_gradients| update_quadrature_points);
  FEValues<dim> fe_values_c (fe_c, qf_cell,
                                   update_values | update_gradients| update_hessians);

  std::vector<Tensor<2, dim> > solution_grads_values (qf_cell.size());
  std::vector<double> solution_c1_values (qf_cell.size());
  std::vector<double> solution_c2_values (qf_cell.size());
  std::vector<double> solution_c3_values (qf_cell.size());

  typename DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  typename DoFHandler<dim>::active_cell_iterator
  cell_c = dof_handler_c.begin_active();
  for (; cell!=endc; ++cell, ++cell_c)
      if (cell->is_locally_owned())
      {
          PointHistory<dim> *lqph =
            reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

          Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
          Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

          Assert(solution_grads_values.size() == n_q_points,
                 ExcInternalError());
          Assert(solution_c1_values.size() == n_q_points,
                 ExcInternalError());

          fe_values.reinit(cell);
          fe_values_c.reinit (cell_c);

          fe_values[u_fe].get_function_gradients(solution,  solution_grads_values);
          fe_values_c.get_function_values(solution_c1,   solution_c1_values);
          fe_values_c.get_function_values(solution_c2,   solution_c2_values);
          fe_values_c.get_function_values(solution_c3,   solution_c3_values);

         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
         {

             lqph[q_point].update_values(solution_grads_values[q_point],
                                         solution_c1_values[q_point],
                                         solution_c2_values[q_point],
                                         solution_c3_values[q_point],
                                         parameters.delta_t,
                                         parameters.L,
                                                fe_values.quadrature_point(q_point),
                                          parameters.a_alpha, parameters.c_alpha, parameters.a_omega, parameters.c_omega);
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
       switch (dim)
         {
         case 1:
           displacement_names.push_back ("displacement");
             break;
         case 2:
             displacement_names.push_back ("x_displacement");
             displacement_names.push_back ("y_displacement");
             break;
         case 3:
             displacement_names.push_back ("x_displacement");
             displacement_names.push_back ("y_displacement");
             displacement_names.push_back ("z_displacement");
             break;
        default:
             Assert (false, ExcNotImplemented());
        }

   data_out.add_data_vector (dof_handler, solution, displacement_names);
   data_out.add_data_vector (dof_handler_c, solution_c0, "c0");
   data_out.add_data_vector (dof_handler_c, solution_c1, "c1");
   data_out.add_data_vector (dof_handler_c, solution_c2, "c2");
   data_out.add_data_vector (dof_handler_c, solution_c3, "c3");

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
std::vector< std::vector< Vector<double> > >
   history_field_stress (dim, std::vector< Vector<double> >(dim)),
   local_history_values_at_qpoints_stress (dim, std::vector< Vector<double> >(dim)),
   local_history_fe_values_stress (dim, std::vector< Vector<double> >(dim));

 for (unsigned int i=0; i<dim; i++)
   for (unsigned int j=0; j<dim; j++)
   {
     history_field_stress[i][j].reinit(history_dof_handler.n_dofs());
     local_history_values_at_qpoints_stress[i][j].reinit(qf_cell.size());
     local_history_fe_values_stress[i][j].reinit(history_fe.dofs_per_cell);
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


 FullMatrix<double> qpoint_to_dof_matrix (history_fe.dofs_per_cell,
                                             qf_cell.size());
 FETools::compute_projection_from_quadrature_points_matrix
           (history_fe,
               qf_cell, qf_cell,
            qpoint_to_dof_matrix);

 typename DoFHandler<dim>::active_cell_iterator
 cell = dof_handler.begin_active(),
 endc = dof_handler.end(),
 dg_cell = history_dof_handler.begin_active();
 for (; cell!=endc; ++cell, ++dg_cell)
   if (cell->is_locally_owned())
   {
     PointHistory<dim> *lqph
            = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
     Assert (lqph >=
                 &quadrature_point_history.front(),
                 ExcInternalError());
     Assert (lqph <
                 &quadrature_point_history.back(),
                 ExcInternalError());
     for (unsigned int i=0; i<dim; i++)
       for (unsigned int j=0; j<dim; j++)
       {
         for (unsigned int q=0; q<qf_cell.size(); ++q)
           {

            local_history_values_at_qpoints_stress[i][j](q) = (lqph[q].get_tau()[i][j])/(lqph[q].get_det_F());
            qpoint_to_dof_matrix.vmult (local_history_fe_values_stress[i][j],
                                        local_history_values_at_qpoints_stress[i][j]);
            dg_cell->set_dof_values (local_history_fe_values_stress[i][j],
                                     history_field_stress[i][j]);

            local_history_values_at_qpoints_drivingforce_c1(q) = lqph[q].get_update_c1();
            qpoint_to_dof_matrix.vmult (local_history_fe_values_drivingforce_c1,
                                        local_history_values_at_qpoints_drivingforce_c1);
            dg_cell->set_dof_values (local_history_fe_values_drivingforce_c1,
                                     history_field_drivingforce_c1);

            local_history_values_at_qpoints_drivingforce_c2(q) = lqph[q].get_update_c2();
            qpoint_to_dof_matrix.vmult (local_history_fe_values_drivingforce_c2,
                                        local_history_values_at_qpoints_drivingforce_c2);
            dg_cell->set_dof_values (local_history_fe_values_drivingforce_c2,
                                     history_field_drivingforce_c2);

            local_history_values_at_qpoints_drivingforce_c3(q) = lqph[q].get_update_c3();
            qpoint_to_dof_matrix.vmult (local_history_fe_values_drivingforce_c3,
                                        local_history_values_at_qpoints_drivingforce_c3);
            dg_cell->set_dof_values (local_history_fe_values_drivingforce_c3,
                                     history_field_drivingforce_c3);
       }
      }
   }

 std::vector<DataComponentInterpretation::DataComponentInterpretation>
               data_component_interpretation2(1, DataComponentInterpretation::component_is_scalar);


 data_out.add_data_vector(history_dof_handler, history_field_stress[0][0], "sigma_11",
                              data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_stress[1][1], "sigma_22",
                              data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_stress[2][2], "sigma_33",
                              data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_stress[0][1], "sigma_12",
                              data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_stress[0][2], "sigma_13",
                              data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_stress[1][2], "sigma_23",
                              data_component_interpretation2);


 data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c1, "dc1",
                                 data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c2, "dc2",
                                    data_component_interpretation2);
 data_out.add_data_vector(history_dof_handler, history_field_drivingforce_c3, "dc3",
                                    data_component_interpretation2);
//////////////////////////
// writing output files
  MappingQEulerian<dim, vectorType > q_mapping(degree, dof_handler, solution);

  Vector<float> subdomain (triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");

  data_out.build_patches(q_mapping, degree);

  const unsigned int cycle = time.get_timestep();

  const std::string filename = ("solution-" +
                                Utilities::int_to_string (cycle, 2) +
                                "." +
                                Utilities::int_to_string
                                (triangulation.locally_owned_subdomain(), 4));
  std::ofstream output ((filename + ".vtu").c_str());
  data_out.write_vtu (output);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
           i<Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
           {
       filenames.push_back ("solution-" +
                             Utilities::int_to_string (cycle, 2) +
                             "." +
                             Utilities::int_to_string (i, 4) +
                             ".vtu");
                           }
       std::ofstream master_output (("solution-" +
                                    Utilities::int_to_string (cycle, 2) +
                                    ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }

 }

template <int dim>
void Solid<dim>::output_resultant_stress()
{

       FEFaceValues<dim> fe_face_values (fe, qf_face,
                                            update_values         | update_quadrature_points  |
                                            update_normal_vectors | update_JxW_values);
          double resultant_force = 0.0;
          double resultant_pseudo_force = 0.0;
          double current_face_area = 0.0;
          double referrence_face_area = 0.0;
          double resultant_E = 0.0;
          typename DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned())
          {

              PointHistory<dim> *lqph =
                                        reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

               for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)

                 if (cell->face(face)->at_boundary()
                    &&
                    (cell->face(face)->boundary_id() == 5))
                  {
                    fe_face_values.reinit (cell, face);
                    for (unsigned int q_point=0; q_point<n_q_points_f; ++q_point)
                      {
                        const Tensor<2, dim> F_inv_tr = lqph[q_point].get_F_inv_tr();
                        const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();
                        const Tensor<2, dim> tau      = lqph[q_point].get_tau();
                        const Tensor<2, dim> E      = lqph[q_point].get_E();
                        const double J = lqph[q_point].get_det_F();

                        const Tensor<1, dim> element_force
                                = (tau*F_inv_tr)* fe_face_values.normal_vector(q_point);
                        const Tensor<1, dim> element_pseudo_force
                                =F_inv*element_force;

                        const double current_element_area = J * ((F_inv_tr*fe_face_values.normal_vector(q_point)).norm());

                        resultant_force += element_force [2]*fe_face_values.JxW(q_point);
                        resultant_pseudo_force += element_pseudo_force [2]*fe_face_values.JxW(q_point);
                          current_face_area += current_element_area* fe_face_values.JxW(q_point);
                          referrence_face_area += fe_face_values.JxW(q_point);
                          resultant_E += E[2][2]*fe_face_values.JxW(q_point);
                      }
                  }
          }


           resultant_cauchy_stress[time.get_timestep()]=-1*Utilities::MPI::sum(resultant_force, mpi_communicator )/
                                                   Utilities::MPI::sum(current_face_area, mpi_communicator );
           resultant_first_piola_stress[time.get_timestep()]=-1*Utilities::MPI::sum(resultant_force, mpi_communicator )/
                                                  Utilities::MPI::sum(referrence_face_area, mpi_communicator );
           resultant_second_piola_stress[time.get_timestep()]=-1*Utilities::MPI::sum(resultant_pseudo_force, mpi_communicator )/
                                                  Utilities::MPI::sum(referrence_face_area, mpi_communicator );
           resultant_lagrangian_strain[time.get_timestep()]=-1*Utilities::MPI::sum(resultant_E, mpi_communicator )/
                                                       Utilities::MPI::sum(referrence_face_area, mpi_communicator );

       if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
         {

           std::ofstream myfile_1;
           std::ofstream myfile_2;
           std::ofstream myfile_3;
           std::ofstream myfile_4;
           std::ofstream myfile_5;
           myfile_1.open ("resultant_cauchy_stress.txt");
           myfile_2.open ("resultant_first_piola_stress.txt");
           myfile_3.open ("resultant_second_piola_stress.txt");
           myfile_4.open ("resultant_lagrangian_strain.txt");
           myfile_5.open ("order_parametre.txt");
           for(unsigned int n=0; n<time.get_timestep(); n++)
           {
             myfile_1<< resultant_cauchy_stress[n]<<std::endl;
             myfile_2<< resultant_first_piola_stress[n]<<std::endl;
             myfile_3<< resultant_second_piola_stress[n]<<std::endl;
              myfile_4<< resultant_lagrangian_strain[n]<<std::endl;
             myfile_5<< order_parameter[n]<<std::endl;
           }
             myfile_1.close();
             myfile_2.close();
             myfile_3.close();
             myfile_4.close();
             myfile_5.close();

      }
   }


//solve
template <int dim>
unsigned int
Solid<dim>::solve ()
{

 vectorType
 completely_distributed_solution (locally_owned_dofs, mpi_communicator);

 SolverControl solver_control (dof_handler.n_dofs(), 1e-6*system_rhs.l2_norm());

 TrilinosWrappers::SolverCG solver(solver_control);


 TrilinosWrappers::PreconditionAMG preconditioner;
 //TrilinosWrappers::PreconditionSSOR preconditioner;
 preconditioner.initialize(tangent_matrix);
 solver.solve (tangent_matrix, completely_distributed_solution, system_rhs,
                     preconditioner);

//      TrilinosWrappers::SolverDirect::AdditionalData data;
//      data.solver_type= "Amesos_Superludist";
//      TrilinosWrappers::SolverDirect solver (solver_control, data);
//      solver.solve (tangent_matrix, completely_distributed_solution, system_rhs);

 constraints.distribute (completely_distributed_solution);

 solution_update = completely_distributed_solution;

 return solver_control.last_step();
}

//Solve phase-field kinetic equation
template <int dim>
void
Solid<dim>::solve_c ( )
{
 assemble_system_c ();
  vectorType   solution_update_c1 (locally_owned_dofs_c, mpi_communicator);
  vectorType   solution_update_c2 (locally_owned_dofs_c, mpi_communicator);
  vectorType   solution_update_c3 (locally_owned_dofs_c, mpi_communicator);

  vectorType   temp_solution_c1 (locally_owned_dofs_c, mpi_communicator);
  vectorType   temp_solution_c2 (locally_owned_dofs_c, mpi_communicator);
  vectorType   temp_solution_c3 (locally_owned_dofs_c, mpi_communicator);

  temp_solution_c1= solution_c1;
  temp_solution_c2= solution_c2;
  temp_solution_c3= solution_c3;


  const double success_tol_c1 = (system_rhs_c1.l2_norm()==0)? 1.e-6 : 1e-6*system_rhs_c1.l2_norm();
  const double success_tol_c2 = (system_rhs_c2.l2_norm()==0)? 1.e-6 : 1e-6*system_rhs_c2.l2_norm();
  const double success_tol_c3 = (system_rhs_c3.l2_norm()==0)? 1.e-6 : 1e-6*system_rhs_c3.l2_norm();

  SolverControl solver_control_c1 (dof_handler_c.n_dofs(), success_tol_c1);
  SolverControl solver_control_c2 (dof_handler_c.n_dofs(), success_tol_c2);
  SolverControl solver_control_c3 (dof_handler_c.n_dofs(), success_tol_c3);

  TrilinosWrappers::SolverCG solver_c1 (solver_control_c1);
  TrilinosWrappers::SolverCG solver_c2 (solver_control_c2);
  TrilinosWrappers::SolverCG solver_c3 (solver_control_c3);

  TrilinosWrappers::PreconditionAMG preconditioner;
  preconditioner.initialize(mass_matrix);

  solver_c1.solve (mass_matrix, solution_update_c1, system_rhs_c1, preconditioner);
  solver_c2.solve (mass_matrix, solution_update_c2, system_rhs_c2, preconditioner);
  solver_c3.solve (mass_matrix, solution_update_c3, system_rhs_c3, preconditioner);

  constraints_c.distribute (solution_update_c1);
  constraints_c.distribute (solution_update_c2);
  constraints_c.distribute (solution_update_c3);

  temp_solution_c1 += solution_update_c1;
  temp_solution_c2 += solution_update_c2;
  temp_solution_c3 += solution_update_c3;

  IndexSet::ElementIterator it=locally_owned_dofs_c.begin()  ;
  unsigned int n_dofs_per_core = locally_owned_dofs_c.n_elements();

  for(unsigned int i = 0; i < n_dofs_per_core; i++)
     {
       if (temp_solution_c1(*it)<  0)
           temp_solution_c1(*it) = 0;
       if (temp_solution_c1(*it)>  1)
           temp_solution_c1(*it) = 1;

       if (temp_solution_c2(*it)<  0)
           temp_solution_c2(*it) = 0;
       if (temp_solution_c2(*it)>  1)
           temp_solution_c2(*it) = 1;

       if (temp_solution_c3(*it)<  0)
           temp_solution_c3(*it) = 0;
       if (temp_solution_c3(*it)>  1)
           temp_solution_c3(*it) = 1;


       it++;


     }

  solution_c1= temp_solution_c1;
  solution_c2= temp_solution_c2;
  solution_c3= temp_solution_c3;

  update_qph_incremental();

  vectorType   temp_solution_c0 (locally_owned_dofs_c, mpi_communicator);
  temp_solution_c0=1;
  temp_solution_c0 -= (temp_solution_c1+temp_solution_c2+temp_solution_c3);
  temp_solution_c0.compress(VectorOperation::add);
  solution_c0= temp_solution_c0;

  pcout << "Solving for Volume fractions: \n"
        << solver_control_c1.last_step()<<"+"
        <<solver_control_c1.last_step()<<"+"
        <<solver_control_c1.last_step()
        << "   CG Solver iterations for C1, C2 and C3."
          << std::endl;


}

// Solve Newton-Raphson iterative alqorithm to solve nonlinear mechanical problem
template <int dim>
void Solid<dim>::solve_nonlinear_timestep()
{
  double initial_rhs_norm = 0.;
  unsigned int newton_iteration = 0;
  unsigned int n_iterations=0;
  vectorType  temp_solution_update(locally_owned_dofs, mpi_communicator);
  vectorType  tmp(locally_owned_dofs, mpi_communicator);
  tmp = solution;
    for (; newton_iteration < 100;   ++newton_iteration)
      {
        make_constraints(newton_iteration);
        assemble_system ();


          if (newton_iteration == 0){
         initial_rhs_norm = system_rhs.l2_norm();
         pcout << " Solving for Displacement:   " << std::endl;
       }
          pcout<<"   rhs_norm : "<<system_rhs.l2_norm();

       // tangent_matrix.print(pcout) ;
       // system_rhs.print(pcout);


         n_iterations = solve ();
         pcout << "    Number of CG iterations: " << n_iterations<< std::endl;

          temp_solution_update = solution_update;


          tmp += temp_solution_update;
          solution = tmp;

          update_qph_incremental();




         if (newton_iteration > 0 && system_rhs.l2_norm() <= 1e-4 * initial_rhs_norm)
          {
           pcout << "CONVERGED! " << std::endl;
           break;
          }
        AssertThrow (newton_iteration < 99,
        ExcMessage("No convergence in nonlinear solver!"));

     }
}


//run
template <int dim>
void Solid<dim>::run()
{
make_grid(); // generates the geometry and mesh
system_setup(); // sets up the system matrices and RHS

// Applying initial condition for volume fraction c
vectorType  tmp_solution_c1(locally_owned_dofs_c, mpi_communicator);
vectorType  tmp_solution_c2(locally_owned_dofs_c, mpi_communicator);
vectorType  tmp_solution_c3(locally_owned_dofs_c, mpi_communicator);
VectorTools::interpolate(dof_handler_c, InitialValues<dim>(1,0), tmp_solution_c1); //initial c
VectorTools::interpolate(dof_handler_c, InitialValues<dim>(2,0), tmp_solution_c2); //initial c
VectorTools::interpolate(dof_handler_c, InitialValues<dim>(3,0), tmp_solution_c3); //initial c
solution_c1= tmp_solution_c1;
solution_c2= tmp_solution_c2;
solution_c3= tmp_solution_c3;

update_qph_incremental();

solve_nonlinear_timestep();
output_results();
output_resultant_stress();

time.increment();

// computed actual time integration to update displacement and c
while (time.current() <= time.end() )
  {

    pcout << std::endl
          << "Time step #" << time.get_timestep() << "; "
          << "advancing to t = " << time.current() << "."
          << std::endl;


//        if (/*time.get_timestep()>3 &&*/ load_step<12 )
//        load_step += 2;

 solve_c();

 solve_nonlinear_timestep();

 if(time.get_timestep()%50 == 0)
  {
   output_results();
  }
 output_resultant_stress();

 time.increment();

  }
}

}
#endif
