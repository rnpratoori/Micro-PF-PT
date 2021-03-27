#ifndef MATERIAL_CONSTITUTIVE_H
#define MATERIAL_CONSTITUTIVE_H

#include "standardtensors.h"


//////////// COMPUTE ELASTIC MODULUS AND STRESSES
namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim>
class Material_Constitutive
{
public:
 Material_Constitutive(const double C_A_11, const double C_A_12, const double C_A_13, const double C_A_33, const double C_A_44,
                        const double C_M_11, const double C_M_12, const double C_M_13, const double C_M_33, const double C_M_44,
                        const double lambda_A_iso,
                                                                       const double mu_A_iso,
                                            const double lambda_M_iso,
                                            const double mu_M_iso,
                                            const double A,
                                            const double delta_psi)
                        // const double a_alpha, const double c_alpha, const double a_omega, const double c_omega)
   :
   det_F(1.0),
   Fe(Tensor<2, dim>()),
   ge(StandardTensors<dim>::I),
   Be(SymmetricTensor<2, dim>()),
   I1(0.0),

   Ge(StandardTensors<dim>::I),
   Ee(Tensor<2, dim>()),
   FeEe(Tensor<2, dim>()),
   EeEe(Tensor<2, dim>()),

   C_A_11(C_A_11),
   C_A_12(C_A_12),
   C_A_13(C_A_13),
   C_A_33(C_A_33),
   C_A_44(C_A_44),
   C_M_11(C_M_11),
   C_M_12(C_M_12),
   C_M_13(C_M_13),
   C_M_33(C_M_33),
   C_M_44(C_M_44),

   lambda_A_iso(lambda_A_iso),
   mu_A_iso(mu_A_iso),
   lambda_M_iso(lambda_M_iso),
   mu_M_iso(mu_M_iso),

   // a_alpha(a_alpha),
   // c_alpha(c_alpha),
   // a_omega(a_omega),
   // c_omega(c_omega),

   C_A(Vector<double> (9)),
   C_M1(Vector<double> (9)),
   C_M2(Vector<double> (9)),
   C_M3(Vector<double> (9)),
   lambda_A(Vector<double> (3)),
   lambda_M1(Vector<double> (3)),
   lambda_M2(Vector<double> (3)),
   lambda_M3(Vector<double> (3)),
   lambda (Vector<double> (3)),
   mu_A(Vector<double> (3)),
   mu_M1(Vector<double> (3)),
   mu_M2(Vector<double> (3)),
   mu_M3(Vector<double> (3)),
   mu (Vector<double> (3)),
   nu_A(Vector<double> (3)),
   nu_M1(Vector<double> (3)),
   nu_M2(Vector<double> (3)),
   nu_M3(Vector<double> (3)),
   nu (Vector<double> (3)),

   A0(A),
   delta_psi0(delta_psi)


 {}

 ~Material_Constitutive()
 {}

 void update_material_data (const Tensor<2, dim> &F, const Tensor<2, dim> &F_e, const double &c_1,const double &c_2,const double &c_3 )
 {
   Fe=F_e;
   det_F = determinant(F);
   ge = symmetrize(Fe*transpose(Fe));
   Ge = symmetrize(transpose(Fe)*Fe);
   Ee = 0.5*(Ge-StandardTensors<dim>::I);
   I1 = trace(Ee);
   FeEe= Fe*Ee;
   EeEe= Ee*Ee;

   c_total=c_1+c_2+c_3;
   c0=1-c_total;
   c1=c_1;
   c2=c_2;
   c3=c_3;

   ///sr=1
   // kd1=0.1753 - 0.0208333 *A0 - 0.0208333 *delta_psi0;
   // kr1=0.149153 + 0.125 *A0 - 0.125 *delta_psi0;
   // kd3=-0.447 + 0.125  *A0 + 0.125 *delta_psi0;
   // kr3=-0.808318 - 1.  *A0 + 1.*delta_psi0;

   kd1=0;
   kr1=0;
   kd3=0;
   kr3=0;

   // Tensor<2, 2*dim> C_Aus;
   // C_Aus[0][0]  = C_Aus[1][1] = C_A_11;
   // C_Aus[2][2]  = C_A_33;
   // C_Aus[3][3]  = C_Aus[4][4] = C_A_44;
   // C_Aus[5][5]  = (C_A_11-C_A_12)/2;
   // C_Aus[0][1]  = C_Aus[1][0] = C_A_12;
   // C_Aus[1][2]  = C_Aus[2][1] = C_Aus[0][2] = C_Aus[2][0] = C_A_13;
   //
   // Tensor<2, 2*dim> C_Mar;
   // C_Mar[0][0]  = C_Mar[1][1] = C_M_11;
   // C_Mar[2][2]  = C_M_33;
   // C_Mar[3][3]  = C_Mar[4][4] = C_M_44;
   // C_Mar[5][5]  = (C_M_11-C_M_12)/2;
   // C_Mar[0][1]  = C_Mar[1][0] = C_M_12;
   // C_Mar[1][2]  = C_Mar[2][1] = C_Mar[0][2] = C_Mar[2][0] = C_M_13;

   // std::cout << C_Aus << '\n';
   // MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

   // std::cout << invert(Rot_mat_4) << '\n';
   // std::cout << C_M2 << '\n';
   // std::cout << C_M3 << '\n';

   // Tensor<2, 2*dim> Rot_mat_4;
   // Rot_mat_4[0][0] = Rot_mat_4[2][2] = 1./4.;
   // Rot_mat_4[0][2] = Rot_mat_4[2][0]  = 3./4.;
   // Rot_mat_4[0][4] = Rot_mat_4[5][3] = sqrt(3.)/2.;
   // Rot_mat_4[1][1] = 1.;
   // Rot_mat_4[2][4] = Rot_mat_4[3][5] = -sqrt(3.)/2.;
   // Rot_mat_4[3][3] = Rot_mat_4[5][5] = 1./2.;
   // Rot_mat_4[0][4] = Rot_mat_4[5][3] = sqrt(3.)/2.;
   // Rot_mat_4[4][0] = -sqrt(3.)/4.;
   // Rot_mat_4[4][2] = sqrt(3.)/4.;
   // Rot_mat_4[4][4] = -1./2.;


// Elstic constants for orthotropic material
   C_A[0] = C_A_11;//C_A_11
   C_A[1] = C_A_11;//C_A_22
   C_A[2] = C_A_33;//C_A_33
   C_A[3] = C_A_44;//C_A_44
   C_A[4] = C_A_44;//C_A_55
   C_A[5] = (C_A_11-C_A_12)/2;//C_A_66
   C_A[6] = C_A_12;//C_A_12
   C_A[7] = C_A_13;//C_A_13
   C_A[8] = C_A_13;//C_A_23



   C_M1[0]  = C_M_11;//C_M1_11
   C_M1[1]  = C_M_11;//C_M1_22
   C_M1[2]  = C_M_33;//C_M1_33
   C_M1[3]  = C_M_44;//C_M1_44
   C_M1[4]  = C_M_44;//C_M1_55
   C_M1[5]  = (C_M_11-C_M_12)/2;//C_M1_66
   C_M1[6]  = C_M_12;//C_M1_12
   C_M1[7]  = C_M_13;//C_M1_13
   C_M1[8]  = C_M_13;//C_M1_23

   // C_M2[0]  = C_M_11;//C_M1_11
   // C_M2[1]  = C_M_11;//C_M1_22
   // C_M2[2]  = C_M_33;//C_M1_33
   // C_M2[3]  = C_M_44;//C_M1_44
   // C_M2[4]  = C_M_44;//C_M1_55
   // C_M2[5]  = (C_M_11-C_M_12)/2;//C_M1_66
   // C_M2[6]  = C_M_12;//C_M1_12
   // C_M2[7]  = C_M_13;//C_M1_13
   // C_M2[8]  = C_M_13;//C_M1_23
   //
   // C_M3[0]  = C_M_11;//C_M1_11
   // C_M3[1]  = C_M_11;//C_M1_22
   // C_M3[2]  = C_M_33;//C_M1_33
   // C_M3[3]  = C_M_44;//C_M1_44
   // C_M3[4]  = C_M_44;//C_M1_55
   // C_M3[5]  = (C_M_11-C_M_12)/2;//C_M1_66
   // C_M3[6]  = C_M_12;//C_M1_12
   // C_M3[7]  = C_M_13;//C_M1_13
   // C_M3[8]  = C_M_13;//C_M1_23

   // C_M2 = Rot_mat_4*C_M1*invert(Rot_mat_4);
   // C_M3 = Rot_mat_4*C_M2*invert(Rot_mat_4);

   C_M2[0]= 137.4;//C_M2_11
   C_M2[1]= 146.7;//C_M2_22
   C_M2[2]= 129.1;//C_M2_33
   C_M2[3]=  35.825;//C_M2_44
   C_M2[4]=  69.5;//C_M2_55
   C_M2[5]=  29.275;//C_M2_66
   C_M2[6]=  70.375;//C_M2_12
   C_M2[7]= 92.75;//C_M2_13
   C_M2[8]=  69.125;//C_M2_23

   C_M3[0]= 137.4;//C_M2_11
   C_M3[1]= 146.7;//C_M2_22
   C_M3[2]= 129.1;//C_M2_33
   C_M3[3]=  35.825;//C_M2_44
   C_M3[4]=  69.5;//C_M2_55
   C_M3[5]=  29.275;//C_M2_66
   C_M3[6]=  70.375;//C_M2_12
   C_M3[7]= 92.75;//C_M2_13
   C_M3[8]=  69.125;//C_M2_23

   // std::cout << Rot_mat_4 << '\n';
   // std::cout << invert(Rot_mat_4) << '\n';
   // std::cout << C_M2 << '\n';
   // std::cout << C_M3 << '\n';
   //
   //
   // MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);

   lambda_A[0]= C_A[0]+C_A[8]+2*C_A[3]-(C_A[6]+C_A[7]+2*C_A[4]+2*C_A[5]);
   lambda_A[1]= C_A[1]+C_A[7]+2*C_A[4]-(C_A[6]+C_A[8]+2*C_A[3]+2*C_A[5]);
   lambda_A[2]= C_A[2]+C_A[6]+2*C_A[5]-(C_A[7]+C_A[8]+2*C_A[3]+2*C_A[4]);

   mu_A[0]= 0.5*(C_A[6]+C_A[7]-C_A[8]);
   mu_A[1]= 0.5*(C_A[6]+C_A[8]-C_A[7]);
   mu_A[2]= 0.5*(C_A[7]+C_A[8]-C_A[6]);

   nu_A[0]= 0.5*(C_A[4]+C_A[5]-C_A[3]);
   nu_A[1]= 0.5*(C_A[3]+C_A[5]-C_A[4]);
   nu_A[2]= 0.5*(C_A[3]+C_A[4]-C_A[5]);

   lambda_M1[0]= C_M1[0]+C_M1[8]+2*C_M1[3]-(C_M1[6]+C_M1[7]+2*C_M1[4]+2*C_M1[5]);
   lambda_M1[1]= C_M1[1]+C_M1[7]+2*C_M1[4]-(C_M1[6]+C_M1[8]+2*C_M1[3]+2*C_M1[5]);
   lambda_M1[2]= C_M1[2]+C_M1[6]+2*C_M1[5]-(C_M1[7]+C_M1[8]+2*C_M1[3]+2*C_M1[4]);

   mu_M1[0]= 0.5*(C_M1[6]+C_M1[7]-C_M1[8]);
   mu_M1[1]= 0.5*(C_M1[6]+C_M1[8]-C_M1[7]);
   mu_M1[2]= 0.5*(C_M1[7]+C_M1[8]-C_M1[6]);

   nu_M1[0]= 0.5*(C_M1[4]+C_M1[5]-C_M1[3]);
   nu_M1[1]= 0.5*(C_M1[3]+C_M1[5]-C_M1[4]);
   nu_M1[2]= 0.5*(C_M1[3]+C_M1[4]-C_M1[5]);

   lambda_M2[0]= C_M2[0]+C_M2[8]+2*C_M2[3]-(C_M2[6]+C_M2[7]+2*C_M2[4]+2*C_M2[5]);
   lambda_M2[1]= C_M2[1]+C_M2[7]+2*C_M2[4]-(C_M2[6]+C_M2[8]+2*C_M2[3]+2*C_M2[5]);
   lambda_M2[2]= C_M2[2]+C_M2[6]+2*C_M2[5]-(C_M2[7]+C_M2[8]+2*C_M2[3]+2*C_M2[4]);

   mu_M2[0]= 0.5*(C_M2[6]+C_M2[7]-C_M2[8]);
   mu_M2[1]= 0.5*(C_M2[6]+C_M2[8]-C_M2[7]);
   mu_M2[2]= 0.5*(C_M2[7]+C_M2[8]-C_M2[6]);

   nu_M2[0]= 0.5*(C_M2[4]+C_M2[5]-C_M2[3]);
   nu_M2[1]= 0.5*(C_M2[3]+C_M2[5]-C_M2[4]);
   nu_M2[2]= 0.5*(C_M2[3]+C_M2[4]-C_M2[5]);

   lambda_M3[0]= C_M3[0]+C_M3[8]+2*C_M3[3]-(C_M3[6]+C_M3[7]+2*C_M3[4]+2*C_M3[5]);
   lambda_M3[1]= C_M3[1]+C_M3[7]+2*C_M3[4]-(C_M3[6]+C_M3[8]+2*C_M3[3]+2*C_M3[5]);
   lambda_M3[2]= C_M3[2]+C_M3[6]+2*C_M3[5]-(C_M3[7]+C_M3[8]+2*C_M3[3]+2*C_M3[4]);

   mu_M3[0]= 0.5*(C_M3[6]+C_M3[7]-C_M3[8]);
   mu_M3[1]= 0.5*(C_M3[6]+C_M3[8]-C_M3[7]);
   mu_M3[2]= 0.5*(C_M3[7]+C_M3[8]-C_M3[6]);

   nu_M3[0]= 0.5*(C_M3[4]+C_M3[5]-C_M3[3]);
   nu_M3[1]= 0.5*(C_M3[3]+C_M3[5]-C_M3[4]);
   nu_M3[2]= 0.5*(C_M3[3]+C_M3[4]-C_M3[5]);

   for (unsigned int n=0; n<3; ++n)
   {
   lambda[n] = lambda_A[n]*c0+lambda_M1[n]*c1+lambda_M2[n]*c2+lambda_M3[n]*c3;
   mu[n] = mu_A[n]*c0+mu_M1[n]*c1+mu_M2[n]*c2+mu_M3[n]*c3;
   nu[n] = nu_A[n]*c0+nu_M1[n]*c1+nu_M2[n]*c2+nu_M3[n]*c3;
   }

// In the case of isotropic material, elstic constants are computed as follows
   lambda_iso = lambda_A_iso+ lambda_M_iso*c_total;
  mu_iso = mu_A_iso+mu_M_iso*c_total;


   Assert(det_F > 0, ExcInternalError());
 }

// Compute the Kirchhoff stress and Jacobians for Orthotropic material
     SymmetricTensor<2, dim> get_tau() const
     {
      SymmetricTensor<2, dim> kirchhoff_stress;
      for (unsigned int n=0; n<dim; ++n)
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)

              kirchhoff_stress[i][j] += lambda[n]*Ee[n][n]*Fe[i][n]*Fe[j][n]+
                                        mu[n]*(I1*Fe[i][n]*Fe[j][n]+Ee[n][n]*ge[i][j])+
                                        2*nu[n]*(Fe[i][n]*FeEe[j][n]+FeEe[i][n]*Fe[j][n]);

     return kirchhoff_stress;
     }
//compute the total Jc
     SymmetricTensor<4, dim> get_Jc() const
     {

       SymmetricTensor<4,dim> elasticityTensor;

          for (unsigned int n=0; n<dim; ++n)
            for (unsigned int i=0; i<dim; ++i)
              for (unsigned int j=0; j<dim; ++j)
                for (unsigned int k=0; k<dim; ++k)
                  for (unsigned int l=0; l<dim; ++l)
                  {
                    elasticityTensor[i][j][k][l] +=
                             lambda[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                 mu[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                 nu[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                        Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
                  }
      return elasticityTensor;
    }
// compute Jc for Austenite and martensitic variants
     SymmetricTensor<4, dim> get_Jc_A() const
       {

         SymmetricTensor<4,dim> elasticityTensor_A;

            for (unsigned int n=0; n<dim; ++n)
              for (unsigned int i=0; i<dim; ++i)
                for (unsigned int j=0; j<dim; ++j)
                  for (unsigned int k=0; k<dim; ++k)
                    for (unsigned int l=0; l<dim; ++l)
                    {
                      elasticityTensor_A[i][j][k][l] +=
                         lambda_A[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                             mu_A[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                             nu_A[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                    Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
                    }
        return elasticityTensor_A;
      }


          SymmetricTensor<4, dim> get_Jc_M1() const
          {

            SymmetricTensor<4,dim> elasticityTensor_M1;

               for (unsigned int n=0; n<dim; ++n)
                 for (unsigned int i=0; i<dim; ++i)
                   for (unsigned int j=0; j<dim; ++j)
                     for (unsigned int k=0; k<dim; ++k)
                       for (unsigned int l=0; l<dim; ++l)
                       {
                         elasticityTensor_M1[i][j][k][l] +=
                            lambda_M1[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                mu_M1[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                nu_M1[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                       Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
                       }
           return elasticityTensor_M1;
         }


             SymmetricTensor<4, dim> get_Jc_M2() const
             {

               SymmetricTensor<4,dim> elasticityTensor_M2;

                  for (unsigned int n=0; n<dim; ++n)
                    for (unsigned int i=0; i<dim; ++i)
                      for (unsigned int j=0; j<dim; ++j)
                        for (unsigned int k=0; k<dim; ++k)
                          for (unsigned int l=0; l<dim; ++l)
                          {
                            elasticityTensor_M2[i][j][k][l] +=
                               lambda_M2[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                   mu_M2[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                   nu_M2[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                          Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
                          }
              return elasticityTensor_M2;
            }


                SymmetricTensor<4, dim> get_Jc_M3() const
                {

                  SymmetricTensor<4,dim> elasticityTensor_M3;

                     for (unsigned int n=0; n<dim; ++n)
                       for (unsigned int i=0; i<dim; ++i)
                         for (unsigned int j=0; j<dim; ++j)
                           for (unsigned int k=0; k<dim; ++k)
                             for (unsigned int l=0; l<dim; ++l)
                             {
                               elasticityTensor_M3[i][j][k][l] +=
                                  lambda_M3[n]*Fe[i][n]*Fe[j][n]*Fe[k][n]*Fe[l][n]+
                                      mu_M3[n]*(Fe[i][n]*Fe[j][n]*ge[k][l]+ ge[i][j]*Fe[k][n]*Fe[l][n])+
                                      nu_M3[n]*(Fe[i][n]*ge[j][k]*Fe[l][n]+ Fe[j][n]*ge[i][k]*Fe[l][n]+
                                             Fe[i][n]*ge[j][l]*Fe[k][n]+ Fe[j][n]*ge[i][l]*Fe[k][n]);
                             }
                 return elasticityTensor_M3;
               }


// compute the driving force excluding the transformational work
 double get_driving_force_noStress () const
 {
  return  det_F*(delta_psi0+A0*(1-2*c_total));
 }

  double get_det_F() const
 {
   return det_F;
 }
// compute the threshhold related terms for the calibration of the instability criteria
 double get_threshold_c1() const
 {
     const double k1=kd1+(kr1-kd1)*c1;
     const double k3=kd3+(kr3-kd3)*c1;

     SymmetricTensor<2, dim>kirchhoff_stress= get_tau();
     const double k_c1=k1*(kirchhoff_stress[0][0]+kirchhoff_stress[1][1])+k3*kirchhoff_stress[2][2];
  return k_c1;
 }

 double get_threshold_c2() const
 {
     const double k1=kd1+(kr1-kd1)*c2;
     const double k3=kd3+(kr3-kd3)*c2;

     SymmetricTensor<2, dim>kirchhoff_stress= get_tau();
     const double k_c2=k1*(kirchhoff_stress[0][0]+kirchhoff_stress[2][2])+k3*kirchhoff_stress[1][1];
  return k_c2;
 }

 double get_threshold_c3() const
 {
     const double k1=kd1+(kr1-kd1)*c3;
     const double k3=kd3+(kr3-kd3)*c3;

     SymmetricTensor<2, dim>kirchhoff_stress= get_tau();
     const double k_c3=k1*(kirchhoff_stress[1][1]+kirchhoff_stress[2][2])+k3*kirchhoff_stress[0][0];
  return k_c3;
 }


protected:
 double det_F;
 // double I1;
 Tensor<2, dim> Fe;
 SymmetricTensor<2, dim> ge;
 SymmetricTensor<2, dim> Be;
 double I1;


 SymmetricTensor<2, dim> Ge;
 Tensor<2, dim> Ee;
 Tensor<2, dim> FeEe;
 Tensor<2, dim> EeEe;

 double C_A_11, C_A_12, C_A_13, C_A_33, C_A_44, C_M_11, C_M_12, C_M_13, C_M_33, C_M_44;
 double lambda_A_iso,mu_A_iso,lambda_M_iso,mu_M_iso,lambda_iso,mu_iso;
 // double a_alpha, c_alpha, a_omega, c_omega;
 // double mu_A_iso,mu_M_iso,mu_iso;

 Vector<double> C_A,C_M1,C_M2,C_M3;
 Vector<double> lambda_A, lambda_M1,lambda_M2,lambda_M3, lambda;
 Vector<double> mu_A, mu_M1,mu_M2,mu_M3,mu;
 Vector<double> nu_A, nu_M1,nu_M2,nu_M3,nu;

 // double lambda_A_iso,lambda_M_iso,lambda_iso;
 // double mu_A_iso,mu_M_iso,mu_iso;
 // double lambda_A_iso,lambda_M_iso,lambda_iso;



 double A0;
 double delta_psi0;
 double ki00;
 double c_total;
 double c0,c1, c2, c3;
 double kd1,kr1,kd3,kr3;



};
}

#endif
