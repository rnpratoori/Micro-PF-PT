# Micro-PF-PT: Micro Phase-Field Phase Transformation Simulation

## Overview

This is a parallel finite element code for simulating phase transformations in materials using a phase-field approach coupled with mechanical deformation. The code is built on the [deal.II](https://www.dealii.org/) finite element library and uses MPI for parallel computation.

## Authors

- Hamed Babaei (2018)
- Raghunandan Pratoori (2020)

## Features

- **Phase-field modeling** of martensitic phase transformations
- **Coupled thermo-mechanical** simulation
- **Parallel computation** using MPI and Trilinos
- **3D finite element** analysis
- **Multiple martensitic variants** (up to 3 variants)
- **Orthotropic material** properties

## Dependencies

- **deal.II** (version 9.2.0 or higher)
- **CMake** (version 2.8.12 or higher)
- **MPI** implementation (e.g., OpenMPI, MPICH)
- **Trilinos** (included with deal.II)

## Building the Code

```bash
# Configure with CMake
cmake .

# Build
make

# The executable 'microPF' will be created
```

## Running Simulations

```bash
# Run with default parameters
./microPF

# The code reads parameters from 'parameters.prm' by default
```

## Code Structure

### Main Components

- **`microPF.cc`**: Main driver program
- **`include/solid.h`**: Main solver class declaration
- **`src/solid.C`**: Main solver class implementation
- **`include/material_constitutive.h`**: Material constitutive model declaration
- **`src/material_constitutive.C`**: Material constitutive model implementation
- **`include/pointhistory.h`**: Quadrature point history declaration
- **`src/pointhistory.C`**: Quadrature point history implementation

### Supporting Files

- **`include/allparameters_str.h`**: Parameter structure aggregation
- **`include/materials_str.h`**: Material properties structure
- **`include/geometry_str.h`**: Geometry parameters
- **`include/fesystem_str.h`**: Finite element system parameters
- **`include/time_str.h`**: Time stepping parameters
- **`include/standardtensors.h`**: Standard tensor definitions
- **`include/boundarydisplacement.h`**: Boundary condition functions
- **`include/initialvalues.h`**: Initial condition functions

## Parameter File

The simulation parameters are specified in `parameters.prm`:

- **Finite element system**: Polynomial degree, quadrature order
- **Geometry**: Mesh refinement, grid scaling
- **Material properties**: Elastic constants for austenite and martensite phases
- **Time stepping**: End time, time step size

## Output Files

The code generates several output files:

- **`solution-*.vtu`**: VTU files for visualization (ParaView compatible)
- **`solution-*.pvtu`**: Parallel VTU files
- **`QuadratureOutputs*.csv`**: Quadrature point data
- **`resultant_*_stress.txt`**: Stress results
- **`resultant_lagrangian_strain.txt`**: Strain results
- **`order_parameter.txt`**: Phase field order parameter evolution

## Code Organization

The code follows a modular structure:

1. **Headers** (`include/`): Class declarations and inline functions
2. **Implementations** (`src/`): Template class implementations in `.C` files
3. **Template inclusion**: Implementation files are included at the end of headers for template instantiation

## License

[Specify license here]

## References

[Add relevant publications or references]

## Contact

For questions or issues, please contact the authors.
