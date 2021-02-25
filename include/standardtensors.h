#ifndef STANDARDTENSORS_H
#define STANDARDTENSORS_H

#include <deal.II/base/symmetric_tensor.h>


//  DEFINE SECOND ORDER IDENTITY, AND TWO FOURTH ORDER IDENTITY TENSORS
namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;
  
template <int dim>
class StandardTensors
{
public:
 static const SymmetricTensor<2, dim> I;
 static const SymmetricTensor<4, dim> IxI;
 static const SymmetricTensor<4, dim> II;
//   static const SymmetricTensor<2, dim> transformation_strain;
};
template <int dim>
const SymmetricTensor<2, dim>
StandardTensors<dim>::I = unit_symmetric_tensor<dim>();
template <int dim>
const SymmetricTensor<4, dim>
StandardTensors<dim>::IxI = outer_product(I, I);
template <int dim>
const SymmetricTensor<4, dim>
StandardTensors<dim>::II = identity_tensor<dim>();
}

#endif
