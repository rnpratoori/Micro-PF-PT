
namespace PhaseField {
template <int dim>
Material_Constitutive<dim>::Material_Constitutive(
    const double C_A_11, const double C_A_12, const double C_A_13,
    const double C_A_33, const double C_A_44, const double C_M_11,
    const double C_M_12, const double C_M_13, const double C_M_33,
    const double C_M_44, const double lambda_A_iso, const double mu_A_iso,
    const double lambda_M_iso, const double mu_M_iso, const double A,
    const double delta_psi)
    : det_F(1.0), Fe(Tensor<2, dim>()), Fe_M2(Tensor<2, dim>()),
      Fe_M3(Tensor<2, dim>()), ge(StandardTensors<dim>::I),
      ge_M2(StandardTensors<dim>::I), ge_M3(StandardTensors<dim>::I),
      Be(SymmetricTensor<2, dim>()), I1(0.0), Ge(StandardTensors<dim>::I),
      Ee(Tensor<2, dim>()), Ee_M2(Tensor<2, dim>()), Ee_M3(Tensor<2, dim>()),
      FeEe(Tensor<2, dim>()), FeEe_M2(Tensor<2, dim>()),
      FeEe_M3(Tensor<2, dim>()), EeEe(Tensor<2, dim>()),
      Rot_mat_2(Tensor<2, dim>()), C_A_11(C_A_11), C_A_12(C_A_12),
      C_A_13(C_A_13), C_A_33(C_A_33), C_A_44(C_A_44), C_M_11(C_M_11),
      C_M_12(C_M_12), C_M_13(C_M_13), C_M_33(C_M_33), C_M_44(C_M_44),
      lambda_A_iso(lambda_A_iso), mu_A_iso(mu_A_iso),
      lambda_M_iso(lambda_M_iso), mu_M_iso(mu_M_iso), C_A(Vector<double>(9)),
      C_M1(Vector<double>(9)), C_M2(Vector<double>(9)), C_M3(Vector<double>(9)),
      lambda_A(Vector<double>(3)), lambda_M1(Vector<double>(3)),
      lambda_M2(Vector<double>(3)), lambda_M3(Vector<double>(3)),
      lambda(Vector<double>(3)), mu_A(Vector<double>(3)),
      mu_M1(Vector<double>(3)), mu_M2(Vector<double>(3)),
      mu_M3(Vector<double>(3)), mu(Vector<double>(3)), nu_A(Vector<double>(3)),
      nu_M1(Vector<double>(3)), nu_M2(Vector<double>(3)),
      nu_M3(Vector<double>(3)), nu(Vector<double>(3)), A0(A),
      delta_psi0(delta_psi) {}

template <int dim> Material_Constitutive<dim>::~Material_Constitutive() {}

template <int dim>
void Material_Constitutive<dim>::update_material_data(const Tensor<2, dim> &F,
                                                      const Tensor<2, dim> &F_e,
                                                      const double &c_1,
                                                      const double &c_2,
                                                      const double &c_3) {
  Fe = F_e;
  det_F = determinant(F);
  ge = symmetrize(Fe * transpose(Fe));
  Ge = symmetrize(transpose(Fe) * Fe);
  Ee = 0.5 * (Ge - StandardTensors<dim>::I);
  I1 = trace(Ee);
  FeEe = Fe * Ee;
  EeEe = Ee * Ee;

  Rot_mat_2[0][0] = Rot_mat_2[1][1] = 1. / 2.;
  Rot_mat_2[0][1] = -sqrt(3.) / 2.;
  Rot_mat_2[1][0] = sqrt(3.) / 2.;
  Rot_mat_2[2][2] = 1.;

  Fe_M2 = Rot_mat_2 * Fe;
  ge_M2 = symmetrize(Fe_M2 * transpose(Fe_M2));
  Ee_M2 =
      0.5 * (symmetrize(transpose(Fe_M2) * Fe_M2) - StandardTensors<dim>::I);
  FeEe_M2 = Fe_M2 * Ee_M2;
  Fe_M3 = Rot_mat_2 * Fe_M2;
  ge_M3 = symmetrize(Fe_M3 * transpose(Fe_M3));
  Ee_M3 =
      0.5 * (symmetrize(transpose(Fe_M3) * Fe_M3) - StandardTensors<dim>::I);
  FeEe_M3 = Fe_M3 * Ee_M3;

  c_total = c_1 + c_2 + c_3;
  c0 = 1 - c_total;
  c1 = c_1;
  c2 = c_2;
  c3 = c_3;

  kd1 = 0.0269978 - A0;

  // Elstic constants for orthotropic material
  C_A[0] = C_A_11;                // C_A_11
  C_A[1] = C_A_11;                // C_A_22
  C_A[2] = C_A_33;                // C_A_33
  C_A[3] = C_A_44;                // C_A_44
  C_A[4] = C_A_44;                // C_A_55
  C_A[5] = (C_A_11 - C_A_12) / 2; // C_A_66
  C_A[6] = C_A_12;                // C_A_12
  C_A[7] = C_A_13;                // C_A_13
  C_A[8] = C_A_13;                // C_A_23

  C_M1[0] = C_M_11;                // C_M1_11
  C_M1[1] = C_M_11;                // C_M1_22
  C_M1[2] = C_M_33;                // C_M1_33
  C_M1[3] = C_M_44;                // C_M1_44
  C_M1[4] = C_M_44;                // C_M1_55
  C_M1[5] = (C_M_11 - C_M_12) / 2; // C_M1_66
  C_M1[6] = C_M_12;                // C_M1_12
  C_M1[7] = C_M_13;                // C_M1_13
  C_M1[8] = C_M_13;                // C_M1_23

  lambda_A[0] = C_A[0] + C_A[8] + 2 * C_A[3] -
                (C_A[6] + C_A[7] + 2 * C_A[4] + 2 * C_A[5]);
  lambda_A[1] = C_A[1] + C_A[7] + 2 * C_A[4] -
                (C_A[6] + C_A[8] + 2 * C_A[3] + 2 * C_A[5]);
  lambda_A[2] = C_A[2] + C_A[6] + 2 * C_A[5] -
                (C_A[7] + C_A[8] + 2 * C_A[3] + 2 * C_A[4]);

  mu_A[0] = 0.5 * (C_A[6] + C_A[7] - C_A[8]);
  mu_A[1] = 0.5 * (C_A[6] + C_A[8] - C_A[7]);
  mu_A[2] = 0.5 * (C_A[7] + C_A[8] - C_A[6]);

  nu_A[0] = 0.5 * (C_A[4] + C_A[5] - C_A[3]);
  nu_A[1] = 0.5 * (C_A[3] + C_A[5] - C_A[4]);
  nu_A[2] = 0.5 * (C_A[3] + C_A[4] - C_A[5]);

  lambda_M1[0] = C_M1[0] + C_M1[8] + 2 * C_M1[3] -
                 (C_M1[6] + C_M1[7] + 2 * C_M1[4] + 2 * C_M1[5]);
  lambda_M1[1] = C_M1[1] + C_M1[7] + 2 * C_M1[4] -
                 (C_M1[6] + C_M1[8] + 2 * C_M1[3] + 2 * C_M1[5]);
  lambda_M1[2] = C_M1[2] + C_M1[6] + 2 * C_M1[5] -
                 (C_M1[7] + C_M1[8] + 2 * C_M1[3] + 2 * C_M1[4]);

  mu_M1[0] = 0.5 * (C_M1[6] + C_M1[7] - C_M1[8]);
  mu_M1[1] = 0.5 * (C_M1[6] + C_M1[8] - C_M1[7]);
  mu_M1[2] = 0.5 * (C_M1[7] + C_M1[8] - C_M1[6]);

  nu_M1[0] = 0.5 * (C_M1[4] + C_M1[5] - C_M1[3]);
  nu_M1[1] = 0.5 * (C_M1[3] + C_M1[5] - C_M1[4]);
  nu_M1[2] = 0.5 * (C_M1[3] + C_M1[4] - C_M1[5]);

  Assert(det_F > 0, ExcInternalError());
}

// Compute the Kirchhoff stress and Jacobians for Orthotropic material
template <int dim>
SymmetricTensor<2, dim> Material_Constitutive<dim>::get_tau() const {
  SymmetricTensor<2, dim> kirchhoff_stress;
  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)

        kirchhoff_stress[i][j] +=
            lambda_A[n] * c0 * Ee[n][n] * Fe[i][n] * Fe[j][n] +
            mu_A[n] * c0 * (I1 * Fe[i][n] * Fe[j][n] + Ee[n][n] * ge[i][j]) +
            2 * nu_A[n] * c0 * (Fe[i][n] * FeEe[j][n] + FeEe[i][n] * Fe[j][n]) +
            lambda_M1[n] * c1 * Ee[n][n] * Fe[i][n] * Fe[j][n] +
            mu_M1[n] * c1 * (I1 * Fe[i][n] * Fe[j][n] + Ee[n][n] * ge[i][j]) +
            2 * nu_M1[n] * c1 *
                (Fe[i][n] * FeEe[j][n] + FeEe[i][n] * Fe[j][n]) +
            lambda_M1[n] * c2 * Ee_M2[n][n] * Fe_M2[i][n] * Fe_M2[j][n] +
            mu_M1[n] * c2 *
                (I1 * Fe_M2[i][n] * Fe_M2[j][n] + Ee_M2[n][n] * ge[i][j]) +
            2 * nu_M1[n] * c2 *
                (Fe_M2[i][n] * FeEe_M2[j][n] + FeEe_M2[i][n] * Fe_M2[j][n]) +
            lambda_M1[n] * c3 * Ee[n][n] * Fe_M2[i][n] * Fe_M2[j][n] +
            mu_M1[n] * c3 *
                (I1 * Fe_M3[i][n] * Fe_M3[j][n] + Ee_M3[n][n] * ge[i][j]) +
            2 * nu_M1[n] * c3 *
                (Fe_M3[i][n] * FeEe_M3[j][n] + FeEe_M3[i][n] * Fe_M3[j][n]);

  return kirchhoff_stress;
}

// compute the total Jc
template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc() const {

  SymmetricTensor<4, dim> elasticityTensor;

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            elasticityTensor[i][j][k][l] +=
                lambda_A[n] * c0 * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_A[n] * c0 *
                    (Fe[i][n] * Fe[j][n] * ge[k][l] +
                     ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_A[n] * c0 *
                    (Fe[i][n] * ge[j][k] * Fe[l][n] +
                     Fe[j][n] * ge[i][k] * Fe[l][n] +
                     Fe[i][n] * ge[j][l] * Fe[k][n] +
                     Fe[j][n] * ge[i][l] * Fe[k][n]) +
                lambda_M1[n] * c1 * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_M1[n] * c1 *
                    (Fe[i][n] * Fe[j][n] * ge[k][l] +
                     ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_M1[n] * c1 *
                    (Fe[i][n] * ge[j][k] * Fe[l][n] +
                     Fe[j][n] * ge[i][k] * Fe[l][n] +
                     Fe[i][n] * ge[j][l] * Fe[k][n] +
                     Fe[j][n] * ge[i][l] * Fe[k][n]) +
                lambda_M1[n] * c2 * Fe_M2[i][n] * Fe_M2[j][n] * Fe_M2[k][n] *
                    Fe_M2[l][n] +
                mu_M1[n] * c2 *
                    (Fe_M2[i][n] * Fe_M2[j][n] * ge_M2[k][l] +
                     ge_M2[i][j] * Fe_M2[k][n] * Fe_M2[l][n]) +
                nu_M1[n] * c2 *
                    (Fe_M2[i][n] * ge_M2[j][k] * Fe_M2[l][n] +
                     Fe_M2[j][n] * ge_M2[i][k] * Fe_M2[l][n] +
                     Fe_M2[i][n] * ge_M2[j][l] * Fe_M2[k][n] +
                     Fe_M2[j][n] * ge_M2[i][l] * Fe_M2[k][n]) +
                lambda_M1[n] * c3 * Fe_M3[i][n] * Fe_M3[j][n] * Fe_M3[k][n] *
                    Fe_M3[l][n] +
                mu_M1[n] * c3 *
                    (Fe_M3[i][n] * Fe_M3[j][n] * ge_M3[k][l] +
                     ge_M3[i][j] * Fe_M3[k][n] * Fe_M3[l][n]) +
                nu_M1[n] * c3 *
                    (Fe_M3[i][n] * ge_M3[j][k] * Fe_M3[l][n] +
                     Fe_M3[j][n] * ge_M3[i][k] * Fe_M3[l][n] +
                     Fe_M3[i][n] * ge_M3[j][l] * Fe_M3[k][n] +
                     Fe_M3[j][n] * ge_M3[i][l] * Fe_M3[k][n]);
          }
  return elasticityTensor;
}

// compute Jc for Austenite and martensitic variants
template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_A() const {

  SymmetricTensor<4, dim> elasticityTensor_A;

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            elasticityTensor_A[i][j][k][l] +=
                lambda_A[n] * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_A[n] * (Fe[i][n] * Fe[j][n] * ge[k][l] +
                           ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_A[n] * (Fe[i][n] * ge[j][k] * Fe[l][n] +
                           Fe[j][n] * ge[i][k] * Fe[l][n] +
                           Fe[i][n] * ge[j][l] * Fe[k][n] +
                           Fe[j][n] * ge[i][l] * Fe[k][n]);
          }
  return elasticityTensor_A;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M1() const {

  SymmetricTensor<4, dim> elasticityTensor_M1;

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            elasticityTensor_M1[i][j][k][l] +=
                lambda_M1[n] * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_M1[n] * (Fe[i][n] * Fe[j][n] * ge[k][l] +
                            ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_M1[n] * (Fe[i][n] * ge[j][k] * Fe[l][n] +
                            Fe[j][n] * ge[i][k] * Fe[l][n] +
                            Fe[i][n] * ge[j][l] * Fe[k][n] +
                            Fe[j][n] * ge[i][l] * Fe[k][n]);
          }
  return elasticityTensor_M1;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M2() const {

  SymmetricTensor<4, dim> elasticityTensor_M2;

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            elasticityTensor_M2[i][j][k][l] +=
                lambda_M1[n] * Fe_M2[i][n] * Fe_M2[j][n] * Fe_M2[k][n] *
                    Fe_M2[l][n] +
                mu_M1[n] * (Fe_M2[i][n] * Fe_M2[j][n] * ge_M2[k][l] +
                            ge_M2[i][j] * Fe_M2[k][n] * Fe_M2[l][n]) +
                nu_M1[n] * (Fe_M2[i][n] * ge_M2[j][k] * Fe_M2[l][n] +
                            Fe_M2[j][n] * ge_M2[i][k] * Fe_M2[l][n] +
                            Fe_M2[i][n] * ge_M2[j][l] * Fe_M2[k][n] +
                            Fe_M2[j][n] * ge_M2[i][l] * Fe_M2[k][n]);
          }
  return elasticityTensor_M2;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M3() const {

  SymmetricTensor<4, dim> elasticityTensor_M3;

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            elasticityTensor_M3[i][j][k][l] +=
                lambda_M1[n] * Fe_M3[i][n] * Fe_M3[j][n] * Fe_M3[k][n] *
                    Fe_M3[l][n] +
                mu_M1[n] * (Fe_M3[i][n] * Fe_M3[j][n] * ge_M3[k][l] +
                            ge_M3[i][j] * Fe_M3[k][n] * Fe_M3[l][n]) +
                nu_M1[n] * (Fe_M3[i][n] * ge_M3[j][k] * Fe_M3[l][n] +
                            Fe_M3[j][n] * ge_M3[i][k] * Fe_M3[l][n] +
                            Fe_M3[i][n] * ge_M3[j][l] * Fe_M3[k][n] +
                            Fe_M3[j][n] * ge_M3[i][l] * Fe_M3[k][n]);
          }
  return elasticityTensor_M3;
}

// compute the driving force excluding the transformational work
template <int dim>
double Material_Constitutive<dim>::get_driving_force_noStress() const {
  return det_F * (delta_psi0 + A0 * (1 - 2 * c_total));
}

template <int dim> double Material_Constitutive<dim>::get_det_F() const {
  return det_F;
}
// compute the threshhold related terms for the calibration of the instability
// criteria
template <int dim> double Material_Constitutive<dim>::get_threshold_c1() const {
  // const double k1=kd1*c1;
  // const double k3=kd3+(kr3-kd3)*c1;

  // SymmetricTensor<2, dim>kirchhoff_stress= get_tau();
  const double k_c1 = kd1;
  return k_c1;
}

template <int dim> double Material_Constitutive<dim>::get_threshold_c2() const {
  const double k1 = kd1 * c2;
  // const double k3=kd3+(kr3-kd3)*c2;

  SymmetricTensor<2, dim> kirchhoff_stress = get_tau();
  const double k_c2 = k1 * (kirchhoff_stress[0][0] + kirchhoff_stress[2][2] +
                            kirchhoff_stress[1][1]);
  return k_c2;
}

template <int dim> double Material_Constitutive<dim>::get_threshold_c3() const {
  const double k1 = kd1 * c3;
  // const double k3=kd3+(kr3-kd3)*c3;

  SymmetricTensor<2, dim> kirchhoff_stress = get_tau();
  const double k_c3 = k1 * (kirchhoff_stress[1][1] + kirchhoff_stress[2][2] +
                            kirchhoff_stress[0][0]);
  return k_c3;
}
} // namespace PhaseField
