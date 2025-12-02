/**
 * @file material_constitutive.h
 * @brief Material constitutive model for phase-transforming materials
 *
 * This class implements the constitutive relations for materials undergoing
 * martensitic phase transformations with orthotropic elastic properties.
 *
 * Key features:
 * - Orthotropic elasticity for austenite and martensite phases
 * - Multiple martensitic variants (up to 3)
 * - Computation of Kirchhoff stress tensor
 * - Computation of material tangent (Jacobian) tensors
 * - Driving forces for phase transformation
 * - Transformation strain calculations
 *
 * Template implementation is in ../src/material_constitutive.C
 */
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
  /**
   * @brief Constructor
   * @param C_A_11 Elastic constant C11 for austenite
   * @param C_A_12 Elastic constant C12 for austenite
   * @param C_A_13 Elastic constant C13 for austenite
   * @param C_A_33 Elastic constant C33 for austenite
   * @param C_A_44 Elastic constant C44 for austenite
   * @param C_M_11 Elastic constant C11 for martensite
   * @param C_M_12 Elastic constant C12 for martensite
   * @param C_M_13 Elastic constant C13 for martensite
   * @param C_M_33 Elastic constant C33 for martensite
   * @param C_M_44 Elastic constant C44 for martensite
   * @param lambda_A_iso Isotropic Lamé parameter λ for austenite
   * @param mu_A_iso Isotropic shear modulus μ for austenite
   * @param lambda_M_iso Isotropic Lamé parameter λ for martensite
   * @param mu_M_iso Isotropic shear modulus μ for martensite
   * @param A Interaction parameter for phase transformation
   * @param delta_psi Chemical free energy difference
   */
  Material_Constitutive(const double C_A_11, const double C_A_12,
                        const double C_A_13, const double C_A_33,
                        const double C_A_44, const double C_M_11,
                        const double C_M_12, const double C_M_13,
                        const double C_M_33, const double C_M_44,
                        const double lambda_A_iso, const double mu_A_iso,
                        const double lambda_M_iso, const double mu_M_iso,
                        const double A, const double delta_psi);

  /** @brief Destructor */
  ~Material_Constitutive();

  /**
   * @brief Update material state based on deformation and phase fields
   * @param F Total deformation gradient
   * @param F_e Elastic deformation gradient
   * @param c_1 Volume fraction of martensite variant 1
   * @param c_2 Volume fraction of martensite variant 2
   * @param c_3 Volume fraction of martensite variant 3
   */
  void update_material_data(const Tensor<2, dim> &F, const Tensor<2, dim> &F_e,
                            const double &c_1, const double &c_2,
                            const double &c_3);

  /**
   * @brief Compute Kirchhoff stress tensor
   * @return Kirchhoff stress τ = J σ
   */
  SymmetricTensor<2, dim> get_tau() const;

  /**
   * @brief Compute total material tangent (Jacobian) tensor
   * @return Fourth-order elasticity tensor Jc
   */
  SymmetricTensor<4, dim> get_Jc() const;

  /** @brief Compute Jacobian for austenite phase */
  SymmetricTensor<4, dim> get_Jc_A() const;

  /** @brief Compute Jacobian for martensite variant 1 */
  SymmetricTensor<4, dim> get_Jc_M1() const;

  /** @brief Compute Jacobian for martensite variant 2 */
  SymmetricTensor<4, dim> get_Jc_M2() const;

  /** @brief Compute Jacobian for martensite variant 3 */
  SymmetricTensor<4, dim> get_Jc_M3() const;

  /**
   * @brief Compute driving force for phase transformation (excluding stress
   * work)
   * @return Chemical driving force
   */
  double get_driving_force_noStress() const;

  /** @brief Get determinant of deformation gradient */
  double get_det_F() const;

  /** @brief Get threshold parameter for variant 1 transformation */
  double get_threshold_c1() const;

  /** @brief Get threshold parameter for variant 2 transformation */
  double get_threshold_c2() const;

  /** @brief Get threshold parameter for variant 3 transformation */
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

  double C_A_11, C_A_12, C_A_13, C_A_33, C_A_44;
  double C_M_11, C_M_12, C_M_13, C_M_33, C_M_44;
  double lambda_A_iso, mu_A_iso;
  double lambda_M_iso, mu_M_iso;

  Vector<double> C_A, C_M1, C_M2, C_M3;
  Vector<double> lambda_A, lambda_M1, lambda_M2, lambda_M3, lambda;
  Vector<double> mu_A, mu_M1, mu_M2, mu_M3, mu;
  Vector<double> nu_A, nu_M1, nu_M2, nu_M3, nu;

  // double lambda_A_iso,lambda_M_iso,lambda_iso;

  double A0;
  double delta_psi0;
  double ki00;
  double c_total;
  double c0, c1, c2, c3;
  double kd1, kr1, kd3, kr3;

  // Cached elasticity tensors for performance (avoid recomputing)
  mutable SymmetricTensor<4, dim> C_A_cached;
  mutable SymmetricTensor<4, dim> C_M1_cached;
  mutable SymmetricTensor<4, dim> C_M2_cached;
  mutable SymmetricTensor<4, dim> C_M3_cached;
  mutable bool tensors_initialized;
};
} // namespace PhaseField

#include "../src/material_constitutive.C"

#endif
