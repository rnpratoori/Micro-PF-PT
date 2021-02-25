#ifndef ALLPARAMETERS_STR_H
#define ALLPARAMETERS_STR_H

#ifndef ALLPARAMETERS_STR_DEF_H
#define ALLPARAMETERS_STR_DEF_H

#include "dealiiheaders.h"
#include "fesystem_str.h"
#include "geometry_str.h"
#include "materials_str.h"
#include "time_str.h"

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

  namespace Parameters
  {
struct AllParameters : public FESystem,
    public Geometry,
    public Materials,
    public Time
  {
    AllParameters(const std::string &input_file);
    static void
    declare_parameters(ParameterHandler &prm);
    void
    parse_parameters(ParameterHandler &prm);
  };
// }
//
// namespace PhaseField
// {
//   using namespace dealii;
//
//   typedef TrilinosWrappers::MPI::Vector vectorType;
//   typedef TrilinosWrappers::SparseMatrix matrixType;
//
//   namespace Parameters
//   {
  // AllParameters::AllParameters(const std::string &input_file)
  // {
  //   ParameterHandler prm;
  //   declare_parameters(prm);
  //   prm.parse_input(input_file);
  //   parse_parameters(prm);
  // }
  // void AllParameters::declare_parameters(ParameterHandler &prm)
  // {
  //   FESystem::declare_parameters(prm);
  //   Geometry::declare_parameters(prm);
  //   Materials::declare_parameters(prm);
  //   Time::declare_parameters(prm);
  // }
  // void AllParameters::parse_parameters(ParameterHandler &prm)
  // {
  //   FESystem::parse_parameters(prm);
  //   Geometry::parse_parameters(prm);
  //   Materials::parse_parameters(prm);
  //   Time::parse_parameters(prm);
  // }
}
}

#endif

#endif
