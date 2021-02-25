#ifndef POINTHISTORY_H
#define POINTHISTORY_H

#include <iostream>
#include <fstream>

#include "dealiiheaders.h"
#include "allparameters_str.h"
#include "standardtensors.h"
#include "material_constitutive.h"


// updates the quadrature point history

namespace PhaseField
{
  using namespace dealii;

  typedef TrilinosWrappers::MPI::Vector vectorType;
  typedef TrilinosWrappers::SparseMatrix matrixType;

template <int dim>
class PointHistory
{
public:
  PointHistory()
    :
    material(NULL),
    F_inv(StandardTensors<dim>::I),
    tau(SymmetricTensor<2, dim>()),
    Jc(SymmetricTensor<4, dim>())
  {}

  virtual ~PointHistory()
  {
    delete material;
    material = NULL;
  }

  void setup_lqp (const Parameters::AllParameters &parameters)
  {
    material = new Material_Constitutive<dim>(parameters.lambdaA,
            parameters.muA,parameters.lambdaM,parameters.muM,parameters.A, parameters.delta_psi);

    update_values(Tensor<2, dim>(), double (),double (),double (), double(), double(), Point<dim>() );
  }

  void update_values (const Tensor<2, dim> &Grad_u_n, const double c1, const double c2, const double c3,
                      const double dt, const double landa, const Point<dim>  q_point )
  {

// Total deformation gradient
    F = (Tensor<2, dim>(StandardTensors<dim>::I) +  Grad_u_n);
// Transformation strain
    Tensor<2, dim> eps_t1;
      eps_t1[0][0] = 0.1753;
     eps_t1[1][1] = 0.1753;
     eps_t1[2][2] = -0.447;

    Tensor<2, dim> eps_t2;
     eps_t2[0][0] = 0.1753;
     eps_t2[1][1] =  -0.447;
     eps_t2[2][2] = 0.1753;

    Tensor<2, dim> eps_t3;
     eps_t3[0][0] = -0.447;
     eps_t3[1][1] = 0.1753;
     eps_t3[2][2] = 0.1753;

    Ft = Tensor<2, dim>(StandardTensors<dim>::I) + eps_t1*c1+eps_t2*c2+eps_t3*c3; // Transformation deformation gradient
    Fe = F * invert(Ft); // Elastic deformation gradient
    material->update_material_data(F, Fe, c1,c2,c3);

    E=0.5*(symmetrize(transpose(F)*F)-StandardTensors<dim>::I); //Total lagrangian strain
    Ee=0.5*(symmetrize(transpose(Fe)*Fe)-StandardTensors<dim>::I); //Elastic lagrangian strain
    F_inv = invert(F);
    F_inv_tr = transpose(F_inv);
    tau = material->get_tau(); // extracting kirchhoff stress
    Jc = material->get_Jc();  // extracting Jacobian
//        Jc_A = material->get_Jc_A();
//        Jc_M1 = material->get_Jc_M1();
//        Jc_M2 = material->get_Jc_M2();
//        Jc_M3 = material->get_Jc_M3();
    driving_force_noStress = material->get_driving_force_noStress(); // extracting driving force with no stress
    k_c1=material->get_threshold_c1();
    k_c2=material->get_threshold_c2();
    k_c3=material->get_threshold_c3();

    const Tensor<2, dim> temp_tensor = F_inv * Tensor<2, dim>(tau);
    const Tensor<2, dim> temp_tensor1 = temp_tensor * Fe;

// driving force from austenite (0) to each martensitic variant (1,2,3) and between the variants.
//The last term can be included if you want to consider change in elasitc property due to PT.
    X10 = scalar_product(temp_tensor1, eps_t1) - driving_force_noStress - k_c1 /*- 0.5*Ee*(Jc_M1-Jc_A)*Ee*/;
    X20 = scalar_product(temp_tensor1, eps_t2) - driving_force_noStress - k_c2 /*- 0.5*Ee*(Jc_M2-Jc_A)*Ee*/;
    X30 = scalar_product(temp_tensor1, eps_t3) - driving_force_noStress - k_c3 /*- 0.5*Ee*(Jc_M3-Jc_A)*Ee*/;

    X21 = scalar_product(temp_tensor1, (eps_t2-eps_t1))/*- 0.5*Ee*(Jc_M2-Jc_M1)*Ee*/;
    X31 = scalar_product(temp_tensor1, (eps_t3-eps_t1))/*- 0.5*Ee*(Jc_M3-Jc_M1)*Ee*/;
    X32 = scalar_product(temp_tensor1, (eps_t3-eps_t2))/*- 0.5*Ee*(Jc_M3-Jc_M2)*Ee*/;

    const double c0=1-c1-c2-c3;

// Implementation of the constraints on the kinetic equation

    if ((X10>0 && c1<1 && c0>0) || (X10<0 && c1>0 && c0<1))
        dc10 = dt * landa * X10;

    if ((X20>0 && c2<1 && c0>0) || (X20<0 && c2>0 && c0<1))
        dc20 = dt * landa * X20;

    if ((X30>0 && c3<1 && c0>0) || (X30<0 && c3>0 && c0<1))
        dc30 = dt * landa * X30;



    if  ((X12>0 && c1<1 && c2>0) || (X12<0 && c1>0 && c2<1))
        dc12 = dt * landa * X12;

    if  ((X13>0 && c1<1 && c3>0) || (X13<0 && c1>0 && c3<1))
        dc13 = dt * landa * X13;

    if  ((X23>0 && c2<1 && c3>0) || (X23<0 && c2>0 && c3<1))
        dc23 = dt * landa * X23;


    if(dc10>0){
        if (c0-dc10<0 || c1+dc10>1)
            dc10=std::min(1-c1,c0);
    }
    else if (dc10<0){
        if (c1-abs(dc10)<0 || c0+abs(dc10)>1)
            dc10=-std::min(c1,1-c0);
    }

    if(dc20>0){
        if (c0-dc20<0 || c2+dc20>1)
            dc20=std::min(1-c2,c0);
    }
    else if (dc20<0){
        if (c2-abs(dc20)<0 || c0+abs(dc20)>1)
            dc20=-std::min(c2,1-c0);
    }

    if(dc30>0){
        if (c0-dc30<0 || c3+dc30>1)
            dc30=std::min(1-c3,c0);
    }
    else if (dc30<0){
        if (c3-abs(dc30)<0 || c0+abs(dc30)>1)
            dc30=-std::min(c3,1-c0);
    }

    if(dc12>0){
        if (c2-dc12<0 || c1+dc12>1)
            dc12=std::min(1-c1,c2);
    }
    else if (dc12<0){
        if (c1-abs(dc12)<0 || c2+abs(dc12)>1)
            dc12=-std::min(c1,1-c2);
    }

    if(dc13>0){
        if (c3-dc13<0 || c1+dc13>1)
            dc13=std::min(1-c1,c3);
    }
    else if (dc13<0){
        if (c1-abs(dc13)<0 || c3+abs(dc13)>1)
            dc13=-std::min(c1,1-c3);
    }

    if(dc23>0){
        if (c3-dc23<0 || c2+dc23>1)
            dc23=std::min(1-c2,c3);
    }
    else if (dc23<0){
        if (c2-abs(dc23)<0 || c3+abs(dc23)>1)
            dc23=-std::min(c2,1-c3);
    }

    const double dci0=dc10+dc20+dc30;

           if(dci0>0){
               if(c0-dci0<0){
                dc10=0;
                   dc20=0;
                   dc30=0;
               }
           }
           else if (dci0<0){
               if(c0+abs(dci0)>1){
                   dc10=0;
                   dc20=0;
                   dc30=0;

               }
           }

    dc21=-dc12;
    dc31=-dc13;
    dc32=-dc23;

    dc1=dc10+dc12+dc13;
    dc2=dc20+dc21+dc23;
    dc3=dc30+dc31+dc32;

// Excluding PT for the top and bottom of the sample to have pure elastic deformation
    if (q_point[2]<0.156 || q_point[2]>2.844){
      dc1=0;
      dc2=0;
      dc3=0;
    }

    Assert(determinant(F_inv) > 0, ExcInternalError());
  }


  const Tensor<2, dim> &get_F() const
  {
    return F;
  }

  double get_det_F() const
  {
    return material->get_det_F();
  }

  const Tensor<2, dim> &get_F_inv() const
  {
    return F_inv;
  }
  const Tensor<2, dim> &get_E() const
  {
    return E;
  }

  const Tensor<2, dim> &get_F_inv_tr() const
  {
    return F_inv_tr;
  }

  const SymmetricTensor<2, dim> &get_tau() const
  {
    return tau;
  }

  const SymmetricTensor<4, dim> &get_Jc() const
  {
    return Jc;
  }

  double get_update_c1() const
  {
    return dc1;
  }
  double get_update_c2() const
  {
    return dc2;
  }
  double get_update_c3() const
  {
    return dc3;
  }

private:
  Material_Constitutive<dim> *material;
  Tensor<2, dim> F;
  Tensor<2, dim> F_inv;
  Tensor<2, dim> F_inv_tr;
  Tensor<2, dim> Ft;
  Tensor<2, dim> E;
  SymmetricTensor<2, dim> Ee;
  SymmetricTensor<2, dim> tau;
  SymmetricTensor<4, dim> Jc, Jc_A,Jc_M1,Jc_M2,Jc_M3;
  Tensor<2, dim> Fe;
  double driving_force_noStress;
  double X10,X12,X13,X20,X21,X23,X30,X31,X32;
  double dc10,dc12,dc13,dc20,dc21,dc23,dc30,dc31,dc32;
  double dc1,dc2,dc3;
  double k_c1,k_c2,k_c3;

};
}
#endif
