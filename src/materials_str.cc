#include "../include/materials_str.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

  namespace Parameters
  {
  void Materials::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Material properties");
    {
      prm.declare_entry("C11 austenite", "144.0", /*  */
                        Patterns::Double(),
                        "C11 austenite");
      prm.declare_entry("C12 austenite", "144.0", /*  */
                        Patterns::Double(),
                        "C12 austenite");
      prm.declare_entry("C13 austenite", "144.0", /*  */
                        Patterns::Double(),
                        "C13 austenite");
      prm.declare_entry("C33 austenite", "144.0", /*  */
                        Patterns::Double(),
                        "C33 austenite");
      prm.declare_entry("C44 austenite", "144.0", /*  */
                        Patterns::Double(),
                        "C44 austenite");
      prm.declare_entry("C11 martensite", "144.0", /*  */
                        Patterns::Double(),
                        "C11 martensite");
      prm.declare_entry("C12 martensite", "144.0", /*  */
                        Patterns::Double(),
                        "C12 martensite");
      prm.declare_entry("C13 martensite", "144.0", /*  */
                        Patterns::Double(),
                        "C13 martensite");
      prm.declare_entry("C33 martensite", "144.0", /*  */
                        Patterns::Double(),
                        "C33 martensite");
      prm.declare_entry("C44 martensite", "144.0", /*  */
                        Patterns::Double(),
                        "C44 martensite");

      prm.declare_entry("lambda austenite", "144.0", /*  */
                        Patterns::Double(),
                        "lambda austenite");
      prm.declare_entry("mu austenite", "74.0",
                        Patterns::Double(0.0),
                        "mu austenite");
      prm.declare_entry("lambda martensite", "379.0", /*  */
                        Patterns::Double(),
                        "lambda martensite");
      prm.declare_entry("mu martensite", "134.0",
                        Patterns::Double(0.0),
                        "mu martensite");
      prm.declare_entry("kinetic coeff", "2.6",
                        Patterns::Double(0.0),
                        "kinetic coeff");
      prm.declare_entry("interaction parameter", "0.028",
                        Patterns::Double(),
                        "interaction parameter");
      prm.declare_entry("thermal jump", "0.028",
                        Patterns::Double(),
                        "thermal jump");
      prm.declare_entry("threshold", "0.028",
                         Patterns::Double(),
                         "threshold");
      prm.declare_entry("a alpha", "0.028",
                        Patterns::Double(),
                        "a alpha");
      prm.declare_entry("c alpha", "0.028",
                        Patterns::Double(),
                        "c alpha");
      prm.declare_entry("a omega", "0.028",
                        Patterns::Double(),
                        "a omega");
      prm.declare_entry("c omega", "0.028",
                        Patterns::Double(),
                        "c omega");
    }
    prm.leave_subsection();
  }

  void Materials::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection("Material properties");
    {
      C_A_11  = prm.get_double("C11 austenite");
      C_A_12  = prm.get_double("C12 austenite");
      C_A_13  = prm.get_double("C13 austenite");
      C_A_33  = prm.get_double("C33 austenite");
      C_A_44  = prm.get_double("C44 austenite");
      C_M_11  = prm.get_double("C11 martensite");
      C_M_12  = prm.get_double("C12 martensite");
      C_M_13  = prm.get_double("C13 martensite");
      C_M_33  = prm.get_double("C33 martensite");
      C_M_44  = prm.get_double("C44 martensite");
      lambdaA = prm.get_double("lambda austenite");
      muA = prm.get_double("mu austenite");
      lambdaM = prm.get_double("lambda martensite");
      muM = prm.get_double("mu martensite");
      L = prm.get_double("kinetic coeff");
      A = prm.get_double("interaction parameter");
      delta_psi = prm.get_double("thermal jump");
      ki0 = prm.get_double("threshold");
      a_alpha = prm.get_double("a alpha");
      c_alpha = prm.get_double("c alpha");
      a_omega = prm.get_double("a omega");
      c_omega = prm.get_double("c omega");
    }
    prm.leave_subsection();
  }
}
}
