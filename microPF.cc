/**
 * @file microPF.cc
 * @brief Main driver program for micro phase-field phase transformation
 * simulation
 *
 * This program simulates phase transformations in materials using a coupled
 * phase-field and mechanical deformation approach. It uses the deal.II finite
 * element library for spatial discretization and MPI for parallel computation.
 *
 * @author Hamed Babaei (2018)
 * @author Raghunandan Pratoori (2020)
 *
 * The simulation couples:
 * - Mechanical equilibrium equations for displacement field
 * - Phase-field kinetic equations for order parameters (martensitic variants)
 * - Orthotropic material constitutive relations
 *
 * The code reads simulation parameters from 'parameters.prm' and outputs
 * results in VTU format for visualization.
 */
#include <fstream>
#include <iostream>
#include <memory>

#include "include/allparameters_str.h"
#include "include/boundarydisplacement.h"
#include "include/dealiiheaders.h"
#include "include/fesystem_str.h"
#include "include/geometry_str.h"
#include "include/initialvalues.h"
#include "include/material_constitutive.h"
#include "include/materials_str.h"
#include "include/pointhistory.h"
#include "include/solid.h"
#include "include/time_str.h"

#include "include/standardtensors.h"
#include "include/timestep.h"

int main(int argc, char *argv[]) {
  try {
    using namespace dealii;
    using namespace PhaseField;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    std::string param_file;
    param_file = "parameters.prm";

    Solid<3> solid_3d(param_file);
    solid_3d.run();
  } catch (std::exception &exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
  } catch (...) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}
