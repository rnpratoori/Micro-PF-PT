#include "../include/restart_str.h"

namespace PhaseField {
namespace Parameters {
void Restart::declare_parameters(dealii::ParameterHandler &prm) {
  prm.enter_subsection("Restart");
  {
    prm.declare_entry("restart", "false", dealii::Patterns::Bool(),
                      "restart a simulation");
  }
  prm.leave_subsection();
}

void Restart::parse_parameters(dealii::ParameterHandler &prm) {
  prm.enter_subsection("Restart");
  {
    restart = prm.get_bool("restart");
  }
  prm.leave_subsection();
}
} // namespace Parameters
} // namespace PhaseField
