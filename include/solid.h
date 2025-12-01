#ifndef SOLID_CORE_H
#define SOLID_CORE_H

#include <iostream>
#include <fstream>
// #include <string>

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
  void    output_quad();
  void    writeQuadratureOutput(unsigned int _currentIncrement);



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

  std::vector<std::vector<double>> outputQuadrature;

};

}
#include "../src/solid.C"
#endif
