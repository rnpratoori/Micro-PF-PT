#include "../include/fesystem_str.h"

namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

namespace Parameters {

void FESystem::declare_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Finite element system");
  {
    prm.declare_entry("Polynomial degree - displacement", "1",
                      Patterns::Integer(0),
                      "Displacement system polynomial order");
    prm.declare_entry("Quadrature order - displacement", "2",
                      Patterns::Integer(0),
                      "Displacement system Gauss quadrature order");
    prm.declare_entry("Polynomial degree - concentration", "1",
                      Patterns::Integer(0),
                      "Concentration system polynomial order");
    prm.declare_entry("Quadrature order - concentration", "2",
                      Patterns::Integer(0),
                      "Concentration system Gauss quadrature order");
  }
  prm.leave_subsection();
}

void FESystem::parse_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Finite element system");
  {
    poly_degree = prm.get_integer("Polynomial degree - displacement");
    quad_order = prm.get_integer("Quadrature order - displacement");
    poly_degree_c = prm.get_integer("Polynomial degree - concentration");
    quad_order_c = prm.get_integer("Quadrature order - concentration");
  }
  prm.leave_subsection();
}
} // namespace Parameters
} // namespace PhaseField
