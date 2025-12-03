#include "../include/materials_str.h"

namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

namespace Parameters {
void Materials::declare_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Material properties");
  {
    prm.declare_entry("C11 austenite", "0.0", Patterns::Double(0.0),
                      "C11 austenite");
    prm.declare_entry("C12 austenite", "0.0", Patterns::Double(0.0),
                      "C12 austenite");
    prm.declare_entry("C13 austenite", "0.0", Patterns::Double(0.0),
                      "C13 austenite");
    prm.declare_entry("C33 austenite", "0.0", Patterns::Double(0.0),
                      "C33 austenite");
    prm.declare_entry("C44 austenite", "0.0", Patterns::Double(0.0),
                      "C44 austenite");

    prm.declare_entry("C11 martensite", "0.0", Patterns::Double(0.0),
                      "C11 martensite");
    prm.declare_entry("C12 martensite", "0.0", Patterns::Double(0.0),
                      "C12 martensite");
    prm.declare_entry("C13 martensite", "0.0", Patterns::Double(0.0),
                      "C13 martensite");
    prm.declare_entry("C33 martensite", "0.0", Patterns::Double(0.0),
                      "C33 martensite");
    prm.declare_entry("C44 martensite", "0.0", Patterns::Double(0.0),
                      "C44 martensite");

    prm.declare_entry("a alpha", "0.0", Patterns::Double(0.0), "a alpha");
    prm.declare_entry("c alpha", "0.0", Patterns::Double(0.0), "c alpha");
    prm.declare_entry("a omega", "0.0", Patterns::Double(0.0), "a omega");
    prm.declare_entry("c omega", "0.0", Patterns::Double(0.0), "c omega");

    prm.declare_entry("kinetic coeff", "0.0", Patterns::Double(0.0),
                      "kinetic coeff");
    prm.declare_entry("interaction parameter", "0.0", Patterns::Double(0.0),
                      "interaction parameter");
    prm.declare_entry("thermal jump", "0.0", Patterns::Double(0.0),
                      "thermal jump");
    prm.declare_entry("limiting threshold", "0.0", Patterns::Double(0.0),
                      "limiting threshold");
  }
  prm.leave_subsection();
}

void Materials::parse_parameters(ParameterHandler &prm) {
  prm.enter_subsection("Material properties");
  {
    C_A_in.clear();
    C_A_in.push_back(prm.get_double("C11 austenite"));
    C_A_in.push_back(prm.get_double("C12 austenite"));
    C_A_in.push_back(prm.get_double("C13 austenite"));
    C_A_in.push_back(prm.get_double("C33 austenite"));
    C_A_in.push_back(prm.get_double("C44 austenite"));

    C_M_in.clear();
    C_M_in.push_back(prm.get_double("C11 martensite"));
    C_M_in.push_back(prm.get_double("C12 martensite"));
    C_M_in.push_back(prm.get_double("C13 martensite"));
    C_M_in.push_back(prm.get_double("C33 martensite"));
    C_M_in.push_back(prm.get_double("C44 martensite"));

    lattice_param.clear();
    lattice_param.push_back(prm.get_double("a alpha"));
    lattice_param.push_back(prm.get_double("c alpha"));
    lattice_param.push_back(prm.get_double("a omega"));
    lattice_param.push_back(prm.get_double("c omega"));

    L = prm.get_double("kinetic coeff");
    A = prm.get_double("interaction parameter");
    delta_psi = prm.get_double("thermal jump");
    k = prm.get_double("limiting threshold");
  }
  prm.leave_subsection();
}
} // namespace Parameters
} // namespace PhaseField
