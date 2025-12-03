
namespace PhaseField {
template <int dim>
void PointHistory<dim>::setup_lqp(const Parameters::AllParameters &parameters) {
  material = new Material_Constitutive<dim>(parameters);

  update_values(Tensor<2, dim>(), double(), double(), double(), double(),
                double(), Point<dim>(), double(), double(), double(), double());
}

template <int dim>
void PointHistory<dim>::update_values(
    const Tensor<2, dim> &Grad_u_n, const double c1, const double c2,
    const double c3, const double dt, const double landa,
    const Point<dim> q_point, const double a_alpha, const double c_alpha,
    const double a_omega, const double c_omega) {

  // Total deformation gradient
  F = (Tensor<2, dim>(StandardTensors<dim>::I) + Grad_u_n);
  // Transformation strain

  Tensor<2, dim> Rot_mat_2;
  Rot_mat_2[0][0] = Rot_mat_2[2][2] = 1. / 2.;
  Rot_mat_2[0][1] = sqrt(3.) / 2.;
  Rot_mat_2[1][0] = -sqrt(3.) / 2.;
  Rot_mat_2[2][2] = 1.;

  Tensor<2, dim> eps_t1;
  eps_t1[0][0] = (c_omega - a_alpha) / a_alpha;
  eps_t1[2][2] = (a_omega - c_alpha) / c_alpha;
  eps_t1[1][1] = (2. * sqrt(3.) * a_omega - 3. * sqrt(3.) * a_alpha) /
                 (3. * sqrt(3.) * a_alpha);

  Tensor<2, dim> eps_t2;
  eps_t2 = Rot_mat_2 * eps_t1 * invert(Rot_mat_2);
  // eps_t2[0][0] = a_alpha;
  // eps_t2[1][1] =  c_alpha;
  // eps_t2[2][2] = a_alpha;

  Tensor<2, dim> eps_t3;
  eps_t3 = Rot_mat_2 * eps_t2 * invert(Rot_mat_2);
  // eps_t3[0][0] = c_alpha;
  // eps_t3[1][1] = a_alpha;
  // eps_t3[2][2] = a_alpha;

  Ft = Tensor<2, dim>(StandardTensors<dim>::I) + eps_t1 * c1 + eps_t2 * c2 +
       eps_t3 * c3;    // Transformation deformation gradient
  Fe = F * invert(Ft); // Elastic deformation gradient
  material->update_material_data(F, Fe, c1, c2, c3);

  E = 0.5 * (symmetrize(transpose(F) * F) -
             StandardTensors<dim>::I); // Total lagrangian strain
  Ee = 0.5 * (symmetrize(transpose(Fe) * Fe) -
              StandardTensors<dim>::I); // Elastic lagrangian strain
  F_inv = invert(F);
  F_inv_tr = transpose(F_inv);
  tau = material->get_tau(); // extracting kirchhoff stress
  Jc = material->get_Jc();   // extracting Jacobian
                             //        Jc_A = material->get_Jc_A();
                             //        Jc_M1 = material->get_Jc_M1();
                             //        Jc_M2 = material->get_Jc_M2();
                             //        Jc_M3 = material->get_Jc_M3();
  driving_force_noStress =
      material->get_driving_force_noStress(); // extracting driving force with
                                              // no stress
  k_c1 = material->get_threshold_c1();
  // k_c2=material->get_threshold_c2();
  // k_c3=material->get_threshold_c3();
  // k_c1=0.0069978;

  const Tensor<2, dim> temp_tensor = F_inv * Tensor<2, dim>(tau);
  const Tensor<2, dim> temp_tensor1 = temp_tensor * Fe;

  // driving force from austenite (0) to each martensitic variant (1,2,3) and
  // between the variants.
  // The last term can be included if you want to consider change in elasitc
  // property due to PT.
  X10 = scalar_product(temp_tensor1, eps_t1) - driving_force_noStress - k_c1 -
        0.5 * Ee * (Jc_M1 - Jc_A) * Ee;
  X20 = scalar_product(temp_tensor1, eps_t2) - driving_force_noStress - k_c1 -
        0.5 * Ee * (Jc_M2 - Jc_A) * Ee;
  X30 = scalar_product(temp_tensor1, eps_t3) - driving_force_noStress - k_c1 -
        0.5 * Ee * (Jc_M3 - Jc_A) * Ee;

  X12 = scalar_product(temp_tensor1,
                       (eps_t1 - eps_t2)) /*- 0.5*Ee*(Jc_M2-Jc_M1)*Ee*/;
  X13 = scalar_product(temp_tensor1,
                       (eps_t1 - eps_t3)) /*- 0.5*Ee*(Jc_M3-Jc_M1)*Ee*/;
  X23 = scalar_product(temp_tensor1,
                       (eps_t2 - eps_t3)) /*- 0.5*Ee*(Jc_M3-Jc_M2)*Ee*/;

  const double c0 = 1 - c1 - c2 - c3;

  // Implementation of the constraints on the kinetic equation

  if ((X10 > 0 && c1 < 1 && c0 > 0) /* || (X10<0 && c1>0 && c0<1)*/) {
    dc10 = dt * landa * X10;
    // std::cout <<"dc10-"<< dc10 << '\n';
    // std::cout <<"X10-"<< X10 << '\n';
    // MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if ((X20 > 0 && c2 < 1 && c0 > 0) /* || (X20<0 && c2>0 && c0<1)*/)
    dc20 = dt * landa * X20;

  if ((X30 > 0 && c3 < 1 && c0 > 0) /* || (X30<0 && c3>0 && c0<1)*/)
    dc30 = dt * landa * X30;

  if ((X12 > 0 && c1 < 1 && c2 > 0) || (X12 < 0 && c1 > 0 && c2 < 1))
    dc12 = dt * landa * X12;

  if ((X13 > 0 && c1 < 1 && c3 > 0) || (X13 < 0 && c1 > 0 && c3 < 1))
    dc13 = dt * landa * X13;

  if ((X23 > 0 && c2 < 1 && c3 > 0) || (X23 < 0 && c2 > 0 && c3 < 1))
    dc23 = dt * landa * X23;

  const double dc10_old = dc10;
  const double dc20_old = dc20;
  const double dc30_old = dc30;
  const double dc12_old = dc12;
  const double dc13_old = dc13;
  const double dc23_old = dc23;

  const double softFactor = 1.0;

  // Constraints on change in concentration
  if ((c0 - dc10 - dc20 - dc30) < -1e-5) {
    dc10 = softFactor * dc10 * c0 /
           (std::abs(dc10_old) + std::abs(dc20_old) + std::abs(dc30_old));
    dc20 = softFactor * dc20 * c0 /
           (std::abs(dc10_old) + std::abs(dc20_old) + std::abs(dc30_old));
    dc30 = softFactor * dc30 * c0 /
           (std::abs(dc10_old) + std::abs(dc20_old) + std::abs(dc30_old));
  }

  if ((c1 + dc10 + dc12 + dc13) < 0) {
    dc10 = softFactor * dc10 * c1 /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
    dc12 = softFactor * dc12 * c1 /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
    dc13 = softFactor * dc13 * c1 /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
  } else if ((c1 + dc10 + dc12 + dc13) > 1) {
    dc10 = softFactor * dc10 * (1 - c1) /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
    dc12 = softFactor * dc12 * (1 - c1) /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
    dc13 = softFactor * dc13 * (1 - c1) /
           (std::abs(dc10_old) + std::abs(dc12_old) + std::abs(dc13_old));
  }

  if ((c2 + dc20 - dc12 + dc23) < 0) {
    dc20 = softFactor * dc20 * c2 /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
    dc12 = softFactor * dc12 * c2 /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
    dc23 = softFactor * dc23 * c2 /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
  } else if ((c2 + dc20 - dc12 + dc23) > 1) {
    dc20 = softFactor * dc20 * (1 - c2) /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
    dc12 = softFactor * dc12 * (1 - c2) /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
    dc23 = softFactor * dc23 * (1 - c2) /
           (std::abs(dc20_old) + std::abs(dc12_old) + std::abs(dc23_old));
  }

  if ((c3 + dc30 - dc13 - dc23) < 0) {
    dc30 = softFactor * dc30 * c3 /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
    dc13 = softFactor * dc13 * c3 /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
    dc23 = softFactor * dc23 * c3 /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
  } else if ((c3 + dc30 - dc13 - dc23) > 1) {
    dc30 = softFactor * dc30 * (1 - c3) /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
    dc13 = softFactor * dc13 * (1 - c3) /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
    dc23 = softFactor * dc23 * (1 - c3) /
           (std::abs(dc30_old) + std::abs(dc13_old) + std::abs(dc23_old));
  }

  dc21 = -dc12;
  dc31 = -dc13;
  dc32 = -dc23;

  dc1 = dc10 + dc12 + dc13;
  dc2 = dc20 + dc21 + dc23;
  dc3 = dc30 + dc31 + dc32;

  // Excluding PT for the top and bottom of the sample to have pure elastic
  // deformation
  if (q_point[2] < 0.156 || q_point[2] > (1.5 - 0.156)) {
    dc1 = 0;
    dc2 = 0;
    dc3 = 0;
  }

  Assert(determinant(F_inv) > 0, ExcInternalError());
}
} // namespace PhaseField
