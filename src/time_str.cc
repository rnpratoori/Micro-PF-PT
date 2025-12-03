#include "../include/time_str.h"

namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

namespace Parameters {

void Time::declare_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Time");
  {
    prm.declare_entry("end time", "50", Patterns::Double(), "end time");
    prm.declare_entry("time step", "0.01", Patterns::Double(), "time step");
  }
  prm.leave_subsection();
}
void Time::parse_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Time");
  {
    end_time = prm.get_double("end time");
    delta_t = prm.get_double("time step");
  }
  prm.leave_subsection();
}
} // namespace Parameters
} // namespace PhaseField
