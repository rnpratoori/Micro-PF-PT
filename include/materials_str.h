#ifndef MATERIALS_STR_H
#define MATERIALS_STR_H

#include "dealiiheaders.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

  namespace Parameters
  {
      struct Materials
          {
            double lambdaA; // austenite phase
            double muA;     // austenite phase
            double lambdaM; // martensite phase
            double muM;     // martensite phase
            double L;        // interface mobility
            double A;       // parameter for interaction energy
            double delta_psi;  //thermal energy jump
            double ki0;  //threshold

            static void
            declare_parameters(ParameterHandler &prm);
            void
            parse_parameters(ParameterHandler &prm);
          };
          // void Materials::declare_parameters(ParameterHandler &prm)
          // {
          //   prm.enter_subsection("Material properties");
          //   {
          //     prm.declare_entry("lambda austenite", "144.0", /*  */
          //                       Patterns::Double(),
          //                       "lambda austenite");
          //     prm.declare_entry("mu austenite", "74.0",
          //                       Patterns::Double(0.0),
          //                       "mu austenite");
          //     prm.declare_entry("lambda martensite", "379.0", /*  */
          //                       Patterns::Double(),
          //                       "lambda martensite");
          //     prm.declare_entry("mu martensite", "134.0",
          //                       Patterns::Double(0.0),
          //                       "mu martensite");
          //     prm.declare_entry("kinetic coeff", "2.6",
          //                       Patterns::Double(0.0),
          //                       "kinetic coeff");
          //     prm.declare_entry("interaction parameter", "0.028",
          //                       Patterns::Double(),
          //                       "interaction parameter");
          //     prm.declare_entry("thermal jump", "0.028",
          //                       Patterns::Double(),
          //                       "thermal jump");
          //     prm.declare_entry("threshold", "0.028",
          //                        Patterns::Double(),
          //                        "threshold");
          //
          //   }
          //   prm.leave_subsection();
          // }
          //
          // void Materials::parse_parameters(ParameterHandler &prm)
          // {
          //   prm.enter_subsection("Material properties");
          //   {
          //     lambdaA = prm.get_double("lambda austenite");
          //     muA = prm.get_double("mu austenite");
          //     lambdaM = prm.get_double("lambda martensite");
          //     muM = prm.get_double("mu martensite");
          //     L = prm.get_double("kinetic coeff");
          //     A = prm.get_double("interaction parameter");
          //     delta_psi = prm.get_double("thermal jump");
          //     ki0 = prm.get_double("threshold");
          //   }
          //   prm.leave_subsection();
          // }
}
}

#endif
