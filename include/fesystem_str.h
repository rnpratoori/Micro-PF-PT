#ifndef FESYSTEM_STR_H
#define FESYSTEM_STR_H

#include "dealiiheaders.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

// INPUT OF PARAMETERS

// #ifndef FESYSTEM
// #define FESYSTEM


namespace Parameters
{
  struct FESystem
  {
    unsigned int poly_degree;
    unsigned int quad_order;
    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
  };

  // void FESystem::declare_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Finite element system");
  //   {
  //     prm.declare_entry("Polynomial degree", "2",
  //                       Patterns::Integer(0),
  //                       "Displacement system polynomial order");
  //     prm.declare_entry("Quadrature order", "3",
  //                       Patterns::Integer(0),
  //                       "Gauss quadrature order");
  //   }
  //   prm.leave_subsection();
  // }
  //
  // void FESystem::parse_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Finite element system");
  //   {
  //     poly_degree = prm.get_integer("Polynomial degree");
  //     quad_order = prm.get_integer("Quadrature order");
  //   }
  //   prm.leave_subsection();
  // }

// #endif
// ////////////////////////////////////////////////////
//

// /////////////////////////////////////////////////
//

//   /////////////////////////////////////////////////
//
// ///////////////////////////////////////////////////////
//
}
}

#endif
