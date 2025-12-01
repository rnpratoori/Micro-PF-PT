/**
 * @file pointhistory.h
 * @brief Quadrature point history for storing material state
 *
 * This class stores and updates the material state at each quadrature point,
 * including:
 * - Deformation gradients (total, elastic, transformation)
 * - Stress tensors (Kirchhoff stress)
 * - Material tangent tensors
 * - Phase field variables and their evolution
 * - Driving forces for phase transformation
 *
 * The class manages the material constitutive model and enforces kinetic
 * constraints on phase transformation.
 *
 * Template implementation is in ../src/pointhistory.C
 */
#ifndef POINTHISTORY_H
#define POINTHISTORY_H

#include <fstream>
#include <iostream>

#include "allparameters_str.h"
#include "dealiiheaders.h"
#include "material_constitutive.h"
#include "standardtensors.h"

// updates the quadrature point history

namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim> class PointHistory {
public:
  PointHistory()
      : material(NULL), F_inv(StandardTensors<dim>::I),
        tau(SymmetricTensor<2, dim>()), Jc(SymmetricTensor<4, dim>()) {}

  virtual ~PointHistory() {
    delete material;
    material = NULL;
  }

  void setup_lqp(const Parameters::AllParameters &parameters);

  void update_values(const Tensor<2, dim> &Grad_u_n, const double c1,
                     const double c2, const double c3, const double dt,
                     const double landa, const Point<dim> q_point,
                     const double a_alpha, const double c_alpha,
                     const double a_omega, const double c_omega);

  const Tensor<2, dim> &get_F() const { return F; }

  double get_det_F() const { return material->get_det_F(); }

  const Tensor<2, dim> &get_Fe() const { return Fe; }

  const Tensor<2, dim> &get_Ft() const { return Ft; }

  const Tensor<2, dim> &get_F_inv() const { return F_inv; }
  const Tensor<2, dim> &get_E() const { return E; }

  const Tensor<2, dim> &get_F_inv_tr() const { return F_inv_tr; }

  const SymmetricTensor<2, dim> &get_tau() const { return tau; }

  const SymmetricTensor<4, dim> &get_Jc() const { return Jc; }

  double get_update_c1() const { return dc1; }
  double get_update_c2() const { return dc2; }
  double get_update_c3() const { return dc3; }

  double get_X10() const { return X10; }
  double get_X20() const { return X20; }
  double get_X30() const { return X30; }
  double get_X12() const { return X12; }
  double get_X13() const { return X13; }
  double get_X21() const { return X21; }
  double get_X23() const { return X23; }
  double get_X31() const { return X31; }
  double get_X32() const { return X32; }

  double get_dc10() const { return dc10; }
  double get_dc20() const { return dc20; }
  double get_dc30() const { return dc30; }
  double get_dc12() const { return dc12; }
  double get_dc13() const { return dc13; }
  double get_dc21() const { return dc21; }
  double get_dc23() const { return dc23; }
  double get_dc31() const { return dc31; }
  double get_dc32() const { return dc32; }

private:
  Material_Constitutive<dim> *material;
  Tensor<2, dim> F;
  Tensor<2, dim> F_inv;
  Tensor<2, dim> F_inv_tr;
  Tensor<2, dim> Ft;
  Tensor<2, dim> E;
  SymmetricTensor<2, dim> Ee;
  SymmetricTensor<2, dim> tau;
  SymmetricTensor<4, dim> Jc, Jc_A, Jc_M1, Jc_M2, Jc_M3;
  Tensor<2, dim> Fe;
  double driving_force_noStress;
  double X10, X12, X13, X20, X21, X23, X30, X31, X32;
  double dc10, dc12, dc13, dc20, dc21, dc23, dc30, dc31, dc32;
  double dc1, dc2, dc3;
  double k_c1, k_c2, k_c3;
};
} // namespace PhaseField

#include "../src/pointhistory.C"

#endif
