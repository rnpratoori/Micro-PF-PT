#ifndef MATERIAL_CONSTITUTIVE_H
#define MATERIAL_CONSTITUTIVE_H

#include "standardtensors.h"

//////////// COMPUTE ELASTIC MODULUS AND STRESSES
namespace PhaseField {
using namespace dealii;

typedef TrilinosWrappers::MPI::Vector vectorType;
typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim> class Material_Constitutive {
public:
  Material_Constitutive(const double C_A_11, const double C_A_12,
                        const double C_A_13, const double C_A_33,
                        const double C_A_44, const double C_M_11,
                        const double C_M_12, const double C_M_13,
                        const double C_M_33, const double C_M_44,
                        const double lambda_A_iso, const double mu_A_iso,
                        const double lambda_M_iso, const double mu_M_iso,
                        const double A, const double delta_psi);

  ~Material_Constitutive();

  void update_material_data(const Tensor<2, dim> &F, const Tensor<2, dim> &F_e,
                            const double &c_1, const double &c_2,
                            const double &c_3);

  // Compute the Kirchhoff stress and Jacobians for Orthotropic material
  SymmetricTensor<2, dim> get_tau() const;
  // compute the total Jc
  SymmetricTensor<4, dim> get_Jc() const;
  // compute Jc for Austenite and martensitic variants
  SymmetricTensor<4, dim> get_Jc_A() const;

  SymmetricTensor<4, dim> get_Jc_M1() const;

  SymmetricTensor<4, dim> get_Jc_M2() const;

  SymmetricTensor<4, dim> get_Jc_M3() const;

  // compute the driving force excluding the transformational work
  double get_driving_force_noStress() const;

  double get_det_F() const;
  // compute the threshhold related terms for the calibration of the instability
  // criteria
  double get_threshold_c1() const;

  double get_threshold_c2() const;

  double get_threshold_c3() const;

protected:
  double det_F;
  // double I1;
  Tensor<2, dim> Fe;
  Tensor<2, dim> Fe_M2;
  Tensor<2, dim> Fe_M3;
  SymmetricTensor<2, dim> ge;
  SymmetricTensor<2, dim> ge_M2;
  SymmetricTensor<2, dim> ge_M3;
  SymmetricTensor<2, dim> Be;
  double I1;

  SymmetricTensor<2, dim> Ge;
  Tensor<2, dim> Ee;
  Tensor<2, dim> Ee_M2;
  Tensor<2, dim> Ee_M3;
  Tensor<2, dim> FeEe;
  Tensor<2, dim> FeEe_M2;
  Tensor<2, dim> FeEe_M3;
  Tensor<2, dim> EeEe;
  Tensor<2, dim> Rot_mat_2;

  double C_A_11, C_A_12, C_A_13, C_A_33, C_A_44, C_M_11, C_M_12, C_M_13, C_M_33,
      C_M_44;
  double lambda_A_iso, mu_A_iso, lambda_M_iso, mu_M_iso, lambda_iso, mu_iso;
  // double mu_A_iso,mu_M_iso,mu_iso;

  Vector<double> C_A, C_M1, C_M2, C_M3;
  Vector<double> lambda_A, lambda_M1, lambda_M2, lambda_M3, lambda;
  Vector<double> mu_A, mu_M1, mu_M2, mu_M3, mu;
  Vector<double> nu_A, nu_M1, nu_M2, nu_M3, nu;

  // double lambda_A_iso,lambda_M_iso,lambda_iso;
  // double mu_A_iso,mu_M_iso,mu_iso;
  // double lambda_A_iso,lambda_M_iso,lambda_iso;

  double A0;
  double delta_psi0;
  double ki00;
  double c_total;
  double c0, c1, c2, c3;
  double kd1, kr1, kd3, kr3;
};
} // namespace PhaseField

#include "../src/material_constitutive.C"

#endif
