#ifndef BOUNDARYDISPLACEMENT_H
#define BOUNDARYDISPLACEMENT_H

#include "dealiiheaders.h"
#include "initialvalues.h"


//Applying displacement on the boundary if periodic condition is absent.
//Displacement on the bpoundary can be alternatively applied through periodic condition as well.

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim>
  class BoundaryDisplacement :  public Function<dim>
  {
  public:
    BoundaryDisplacement (const double direction, const int timestep);
    virtual
    void
    vector_value (const Point<dim> &p,
                  Vector<double>   &values) const;
    virtual
    void
    vector_value_list (const std::vector<Point<dim> > &points,
                       std::vector<Vector<double> >   &value_list) const;
  private:
    const double direction;
    const int timestep;
  };

template <int dim>
BoundaryDisplacement<dim>::BoundaryDisplacement (const double direction, const int timestep)
  :
  Function<dim> (dim),
  direction(direction),
  timestep(timestep)
{}

template <int dim>
inline
void
BoundaryDisplacement<dim>::vector_value (const Point<dim> &/*p*/,
                              Vector<double>   &values) const
{
  Assert (values.size() == dim,
          ExcDimensionMismatch (values.size(), dim));

  values = 0.0;
  //values(direction)=-5e-4;

  if(timestep<70)
      values(direction)=-1.5e-3;
  else if(timestep%10==0)
      values(direction)=-1.5e-3;
  else
      values(direction)=0.0;
}

template <int dim>
void
BoundaryDisplacement<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                   std::vector<Vector<double> >   &value_list) const
{
  const unsigned int n_ooints = points.size();
  Assert (value_list.size() == n_ooints,
          ExcDimensionMismatch (value_list.size(), n_ooints));
  for (unsigned int p=0; p<n_ooints; ++p)
      BoundaryDisplacement<dim>::vector_value (points[p],
                                  value_list[p]);
}

}
#endif
