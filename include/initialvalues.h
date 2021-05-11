#ifndef INITIALVALUES_H
#define INITIALVALUES_H

#include "dealiiheaders.h"
#include "initialvalues.h"


// defines the initial condition for the order parameter

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim>
class InitialValues : public Function<dim>
{
public:
 InitialValues (const int &variant,const int &time_step)
  :
   Function<dim>(),
     variant (variant),
   time_step (time_step)

    {}
 virtual double value(const Point<dim>   &p,
                       const unsigned int  /*component = 0*/) const;
private:
 const int variant ;
 const int time_step ;

};

template <int dim>
double InitialValues<dim>::value (const Point<dim>  &p,
                           const unsigned int /*component*/) const
{

//               if (pow(p[0]-35,2)+pow(p[1]-10,2)<4)
//                {
//                   if (variant==1)
//                       return 0.01;
//                   else
//                       return 0;
//                }
//               else if (time_step==0)
        return 0.0;

}

}
#endif
