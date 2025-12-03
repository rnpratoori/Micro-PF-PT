#ifndef RESTART_STR_H
#define RESTART_STR_H

#include <deal.II/base/parameter_handler.h>

namespace PhaseField {
namespace Parameters {
/////////////////// Class to declare restart conditions
struct Restart {
  bool restart;

  static void declare_parameters(dealii::ParameterHandler &prm);

  void parse_parameters(dealii::ParameterHandler &prm);
};
} // namespace Parameters
} // namespace PhaseField

#endif // RESTART_STR_H
