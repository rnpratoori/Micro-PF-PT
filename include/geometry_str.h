#ifndef GEOMETRY_STR_H
#define GEOMETRY_STR_H

#include "dealiiheaders.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

  namespace Parameters
  {
  struct Geometry
  {
    unsigned int refinement;
    double       scale;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
  // void Geometry::declare_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Geometry");
  //   {
  //     prm.declare_entry("Global refinement", "4",
  //                        Patterns::Integer(0),
  //                        "Global refinement level");
  //     prm.declare_entry("Grid scale", "1",
  //                       Patterns::Double(0.0),
  //                       "Global grid scaling factor");
  //   }
  //   prm.leave_subsection();
  // }
  // void Geometry::parse_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Geometry");
  //   {
  //     refinement = prm.get_integer("Global refinement");
  //     scale = prm.get_double("Grid scale");
  //   }
  //   prm.leave_subsection();
  // }
}
}

#endif
