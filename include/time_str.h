#ifndef TIME_STR_H
#define TIME_STR_H

#include "dealiiheaders.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

  namespace Parameters
  {

  struct Time
  {
    double delta_t;
    double end_time;
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
  // void Time::declare_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Time");
  //   {
  //     prm.declare_entry("End time", "50",
  //                       Patterns::Double(),
  //                       "End time");
  //     prm.declare_entry("Time step size", "0.01",
  //                       Patterns::Double(),
  //                       "Time step size");
  //   }
  //   prm.leave_subsection();
  // }
  // void Time::parse_parameters(ParameterHandler &prm)
  // {
  //   prm.enter_subsection("Time");
  //   {
  //     end_time = prm.get_double("End time");
  //     delta_t = prm.get_double("Time step size");
  //   }
  //   prm.leave_subsection();
  // }
}
}

#endif
