#include "../include/fesystem_str.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

namespace Parameters
{

void FESystem::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Finite element system");
  {
    prm.declare_entry("Polynomial degree", "2",
                      Patterns::Integer(0),
                      "Displacement system polynomial order");
    prm.declare_entry("Quadrature order", "3",
                      Patterns::Integer(0),
                      "Gauss quadrature order");
  }
  prm.leave_subsection();
}

void FESystem::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("Finite element system");
  {
    poly_degree = prm.get_integer("Polynomial degree");
    quad_order = prm.get_integer("Quadrature order");
  }
  prm.leave_subsection();
}
}
}
