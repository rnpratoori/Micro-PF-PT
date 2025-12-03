
namespace PhaseField {
template <int dim>
Material_Constitutive<dim>::Material_Constitutive(
    const Parameters::Materials &parameters)
    : tensors_initialized(false) {

  // Initialize C_A and C_M from parameters
  C_A.reinit(parameters.C_A_in.size());
  for (unsigned int i = 0; i < parameters.C_A_in.size(); ++i)
    C_A[i] = parameters.C_A_in[i];

  C_M1.reinit(parameters.C_M_in.size());
  for (unsigned int i = 0; i < parameters.C_M_in.size(); ++i)
    C_M1[i] = parameters.C_M_in[i];

  // Initialize Material_Constitutive members from parameters
  this->A0 = parameters.A;
  this->delta_psi0 = parameters.delta_psi;
  this->ki00 = parameters.k;

  // Initialize lambda/mu/nu vectors
  lambda_A.reinit(3);
  mu_A.reinit(3);
  nu_A.reinit(3);
  lambda_M1.reinit(3);
  mu_M1.reinit(3);
  nu_M1.reinit(3);
  lambda_M2.reinit(3);
  mu_M2.reinit(3);
  nu_M2.reinit(3);
  lambda_M3.reinit(3);
  mu_M3.reinit(3);
  nu_M3.reinit(3);
}

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

  kd1 = 0.0269978 - A0;

  // C_A and C_M1 are initialized in constructor from parameters (size 5).
  // We need to map them to 9-element orthotropic representation.
  // CRITICAL FIX: reinit(9) clears the vector, so we must extract BEFORE
  // resizing.

  double C11_A, C12_A, C13_A, C33_A, C44_A;

  if (C_A.size() == 5) {
    C11_A = C_A[0];
    C12_A = C_A[1];
    C13_A = C_A[2];
    C33_A = C_A[3];
    C44_A = C_A[4];

    C_A.reinit(9);
    C_A[0] = C11_A;               // C_A_11
    C_A[1] = C11_A;               // C_A_22
    C_A[2] = C33_A;               // C_A_33
    C_A[3] = C44_A;               // C_A_44
    C_A[4] = C44_A;               // C_A_55
    C_A[5] = (C11_A - C12_A) / 2; // C_A_66
    C_A[6] = C12_A;               // C_A_12
    C_A[7] = C13_A;               // C_A_13
    C_A[8] = C13_A;               // C_A_23
  } else if (C_A.size() == 9) {
    // Already mapped
    C11_A = C_A[0];
    C12_A = C_A[6];
    C13_A = C_A[7];
    C33_A = C_A[2];
    C44_A = C_A[3];
  } else {
    // Fallback or error
    C11_A = C12_A = C13_A = C33_A = C44_A = 0.0;
    if (C_A.size() != 9)
      C_A.reinit(9);
  }

  double C11_M, C12_M, C13_M, C33_M, C44_M;

  if (C_M1.size() == 5) {
    C11_M = C_M1[0];
    C12_M = C_M1[1];
    C13_M = C_M1[2];
    C33_M = C_M1[3];
    C44_M = C_M1[4];

    C_M1.reinit(9);
    C_M1[0] = C11_M;               // C_M1_11
    C_M1[1] = C11_M;               // C_M1_22
    C_M1[2] = C33_M;               // C_M1_33
    C_M1[3] = C44_M;               // C_M1_44
    C_M1[4] = C44_M;               // C_M1_55
    C_M1[5] = (C11_M - C12_M) / 2; // C_M1_66
    C_M1[6] = C12_M;               // C_M1_12
    C_M1[7] = C13_M;               // C_M1_13
    C_M1[8] = C13_M;               // C_M1_23
  } else if (C_M1.size() == 9) {
    // Already mapped
    C11_M = C_M1[0];
    C12_M = C_M1[6];
    C13_M = C_M1[7];
    C33_M = C_M1[2];
    C44_M = C_M1[3];
  } else {
    C11_M = C12_M = C13_M = C33_M = C44_M = 0.0;
    if (C_M1.size() != 9)
      C_M1.reinit(9);
  }

  kd1 = 0.0269978 - A0;

  // Compute elastic constants
  // C_A and C_M are members initialized in constructor.
  // We need to ensure they have correct size if not already.
  // Constructor initializes them from parameters.

  // We assume C_A has 5 elements from input, but we need 9 for calculation?
  // If input only gives 5, we need to map them.
  // C_A[0] = C11, C_A[1] = C12...
  // Let's check how they are stored in C_A vector.
  // In Materials::parse_parameters, we pushed back 5 values.
  // So C_A has size 5.
  // But the code below uses indices up to 8!
  // C_A[8] = C_A_23 which is same as C_A_13.

  // We need to map the 5 input values to the 9 values expected by the logic.
  // Or update the logic to use the 5 input values.

  // Let's create a local helper or just map them.
  // C11 = C_A[0]
  // C12 = C_A[1]
  // C13 = C_A[2]
  // C33 = C_A[3]
  // C44 = C_A[4]
  // Already declared above at lines 86-91

  // Let's populate the members lambda_A, mu_A, nu_A.
  // They are Vector<double> of size 3.
  // We need to ensure they are sized.
  if (lambda_A.size() != 3)
    lambda_A.reinit(3);
  if (mu_A.size() != 3)
    mu_A.reinit(3);
  if (nu_A.size() != 3)
    nu_A.reinit(3);

  if (lambda_M1.size() != 3)
    lambda_M1.reinit(3);
  if (mu_M1.size() != 3)
    mu_M1.reinit(3);
  if (nu_M1.size() != 3)
    nu_M1.reinit(3);

  // Austenite
  // 0: 11, 1: 22, 2: 33
  // C_A indices in original code:
  // 0:11, 1:22, 2:33, 3:44, 4:55, 5:66, 6:12, 7:13, 8:23

  // Mapping from 5 inputs:
  // C11 -> 0, 1
  // C33 -> 2
  // C44 -> 3, 4
  // (C11-C12)/2 -> 5
  // C12 -> 6
  // C13 -> 7, 8

  double c_a_0 = C11_A;
  double c_a_1 = C11_A;
  double c_a_2 = C33_A;
  double c_a_3 = C44_A;
  double c_a_4 = C44_A;
  double c_a_5 = (C11_A - C12_A) / 2.0;
  double c_a_6 = C12_A;
  double c_a_7 = C13_A;
  double c_a_8 = C13_A;

  lambda_A[0] =
      c_a_0 + c_a_8 + 2 * c_a_3 - (c_a_6 + c_a_7 + 2 * c_a_4 + 2 * c_a_5);
  lambda_A[1] =
      c_a_1 + c_a_7 + 2 * c_a_4 - (c_a_6 + c_a_8 + 2 * c_a_3 + 2 * c_a_5);
  lambda_A[2] =
      c_a_2 + c_a_6 + 2 * c_a_5 - (c_a_7 + c_a_8 + 2 * c_a_3 + 2 * c_a_4);

  mu_A[0] = 0.5 * (c_a_6 + c_a_7 - c_a_8);
  mu_A[1] = 0.5 * (c_a_6 + c_a_8 - c_a_7);
  mu_A[2] = 0.5 * (c_a_7 + c_a_8 - c_a_6);

  nu_A[0] = 0.5 * (c_a_4 + c_a_5 - c_a_3);
  nu_A[1] = 0.5 * (c_a_3 + c_a_5 - c_a_4);
  nu_A[2] = 0.5 * (c_a_3 + c_a_4 - c_a_5);

  // Martensite constants already declared above (lines 104-108)

  double c_m_0 = C11_M;
  double c_m_1 = C11_M;
  double c_m_2 = C33_M;
  double c_m_3 = C44_M;
  double c_m_4 = C44_M;
  double c_m_5 = (C11_M - C12_M) / 2.0;
  double c_m_6 = C12_M;
  double c_m_7 = C13_M;
  double c_m_8 = C13_M;

  lambda_M1[0] =
      c_m_0 + c_m_8 + 2 * c_m_3 - (c_m_6 + c_m_7 + 2 * c_m_4 + 2 * c_m_5);
  lambda_M1[1] =
      c_m_1 + c_m_7 + 2 * c_m_4 - (c_m_6 + c_m_8 + 2 * c_m_3 + 2 * c_m_5);
  lambda_M1[2] =
      c_m_2 + c_m_6 + 2 * c_m_5 - (c_m_7 + c_m_8 + 2 * c_m_3 + 2 * c_m_4);

  mu_M1[0] = 0.5 * (c_m_6 + c_m_7 - c_m_8);
  mu_M1[1] = 0.5 * (c_m_6 + c_m_8 - c_m_7);
  mu_M1[2] = 0.5 * (c_m_7 + c_m_8 - c_m_6);

  nu_M1[0] = 0.5 * (c_m_4 + c_m_5 - c_m_3);
  nu_M1[1] = 0.5 * (c_m_3 + c_m_5 - c_m_4);
  nu_M1[2] = 0.5 * (c_m_3 + c_m_4 - c_m_5);

  // Populate M2 and M3 if needed (assuming isotropic-like or rotated
  // properties) For now, just copy M1 to M2 and M3 if they are not set? The
  // full version logic for M2/M3 is complex if they are rotated. But here we
  // compute lambda/mu/nu which are scalar-like components in Voigt notation?
  // No, they are vectors of 3 components.

  // Let's assume M2 and M3 have same properties as M1 for now, as in original
  // code original code used lambda_M_iso etc.

  if (lambda_M2.size() != 3)
    lambda_M2.reinit(3);
  if (mu_M2.size() != 3)
    mu_M2.reinit(3);
  if (nu_M2.size() != 3)
    nu_M2.reinit(3);

  if (lambda_M3.size() != 3)
    lambda_M3.reinit(3);
  if (mu_M3.size() != 3)
    mu_M3.reinit(3);
  if (nu_M3.size() != 3)
    nu_M3.reinit(3);

  lambda_M2 = lambda_M1;
  mu_M2 = mu_M1;
  nu_M2 = nu_M1;

  lambda_M3 = lambda_M1;
  mu_M3 = mu_M1;
  nu_M3 = nu_M1;

  // Invalidate cached elasticity tensors since Fe has changed
  tensors_initialized = false;

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

  // Return cached value if already computed
  if (tensors_initialized) {
    return C_A_cached;
  }

  // Compute all four tensors once and cache them
  C_A_cached = SymmetricTensor<4, dim>();
  C_M1_cached = SymmetricTensor<4, dim>();
  C_M2_cached = SymmetricTensor<4, dim>();
  C_M3_cached = SymmetricTensor<4, dim>();

  for (unsigned int n = 0; n < dim; ++n)
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        for (unsigned int k = 0; k < dim; ++k)
          for (unsigned int l = 0; l < dim; ++l) {
            C_A_cached[i][j][k][l] +=
                lambda_A[n] * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_A[n] * (Fe[i][n] * Fe[j][n] * ge[k][l] +
                           ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_A[n] * (Fe[i][n] * ge[j][k] * Fe[l][n] +
                           Fe[j][n] * ge[i][k] * Fe[l][n] +
                           Fe[i][n] * ge[j][l] * Fe[k][n] +
                           Fe[j][n] * ge[i][l] * Fe[k][n]);

            C_M1_cached[i][j][k][l] +=
                lambda_M1[n] * Fe[i][n] * Fe[j][n] * Fe[k][n] * Fe[l][n] +
                mu_M1[n] * (Fe[i][n] * Fe[j][n] * ge[k][l] +
                            ge[i][j] * Fe[k][n] * Fe[l][n]) +
                nu_M1[n] * (Fe[i][n] * ge[j][k] * Fe[l][n] +
                            Fe[j][n] * ge[i][k] * Fe[l][n] +
                            Fe[i][n] * ge[j][l] * Fe[k][n] +
                            Fe[j][n] * ge[i][l] * Fe[k][n]);

            C_M2_cached[i][j][k][l] +=
                lambda_M2[n] * Fe_M2[i][n] * Fe_M2[j][n] * Fe_M2[k][n] *
                    Fe_M2[l][n] +
                mu_M2[n] * (Fe_M2[i][n] * Fe_M2[j][n] * ge_M2[k][l] +
                            ge_M2[i][j] * Fe_M2[k][n] * Fe_M2[l][n]) +
                nu_M2[n] * (Fe_M2[i][n] * ge_M2[j][k] * Fe_M2[l][n] +
                            Fe_M2[j][n] * ge_M2[i][k] * Fe_M2[l][n] +
                            Fe_M2[i][n] * ge_M2[j][l] * Fe_M2[k][n] +
                            Fe_M2[j][n] * ge_M2[i][l] * Fe_M2[k][n]);

            C_M3_cached[i][j][k][l] +=
                lambda_M3[n] * Fe_M3[i][n] * Fe_M3[j][n] * Fe_M3[k][n] *
                    Fe_M3[l][n] +
                mu_M3[n] * (Fe_M3[i][n] * Fe_M3[j][n] * ge_M3[k][l] +
                            ge_M3[i][j] * Fe_M3[k][n] * Fe_M3[l][n]) +
                nu_M3[n] * (Fe_M3[i][n] * ge_M3[j][k] * Fe_M3[l][n] +
                            Fe_M3[j][n] * ge_M3[i][k] * Fe_M3[l][n] +
                            Fe_M3[i][n] * ge_M3[j][l] * Fe_M3[k][n] +
                            Fe_M3[j][n] * ge_M3[i][l] * Fe_M3[k][n]);
          }

  tensors_initialized = true;
  return C_A_cached;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M1() const {
  // Ensure cache is initialized by calling get_Jc_A if needed
  if (!tensors_initialized) {
    get_Jc_A(); // This will compute and cache all four tensors
  }
  return C_M1_cached;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M2() const {
  // Ensure cache is initialized
  if (!tensors_initialized) {
    get_Jc_A(); // This will compute and cache all four tensors
  }
  return C_M2_cached;
}

template <int dim>
SymmetricTensor<4, dim> Material_Constitutive<dim>::get_Jc_M3() const {
  // Ensure cache is initialized
  if (!tensors_initialized) {
    get_Jc_A(); // This will compute and cache all four tensors
  }
  return C_M3_cached;
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
