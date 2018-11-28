//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADTemperature.h"

registerADMooseObject("NavierStokesApp", INSADTemperature);

defineADValidParams(INSADTemperature,
                    INSADBase,
                    params.addClassDescription(
                        "This class computes the residual and Jacobian contributions for the "
                        "incompressible Navier-Stokes temperature (energy) equation.");
                    // Coupled variables
                    params.addRequiredCoupledVar("u", "x-velocity");
                    params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
                    params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

                    // Optional parameters
                    params.addParam<MaterialPropertyName>("rho_name", "rho", "density name");
                    params.addParam<MaterialPropertyName>("k_name",
                                                          "k",
                                                          "thermal conductivity name");
                    params.addParam<MaterialPropertyName>("cp_name", "cp", "specific heat name");
                    params.addParam<MaterialPropertyName>("grad_k", 0, "The gradient of k");
                    params.addParam<bool>("supg", false, "Whether to perform SUPG stabilization"););

template <ComputeStage compute_stage>
INSADTemperature<compute_stage>::INSADTemperature(const InputParameters & parameters)
  : INSADBase<compute_stage>(parameters),

    // Material Properties
    _k(adGetADMaterialProperty<Real>("k_name")),
    _cp(adGetADMaterialProperty<Real>("cp_name")),
    _grad_k(adGetADMaterialProperty<RealVectorValue>("grad_k")),
    _second_u(_var.template adSecondSln<compute_stage>()),
    _u_dot(_var.template adUDot<compute_stage>()),
    _supg(adGetParam<bool>("supg"))
{
}

template <ComputeStage compute_stage>
ADResidual
INSADTemperature<compute_stage>::computeQpResidual()
{
  INSVectorValue<compute_stage> U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  auto strong_convective_part = _rho[_qp] * _cp[_qp] * U * _grad_u[_qp];

  auto residual =
      strong_convective_part * _test[_i][_qp] + _k[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  if (_supg)
  {
    residual += strong_convective_part * this->tau() * U * _grad_test[_i][_qp];
    // residual += -(_grad_k[_qp] * _grad_u[_qp] + _k[_qp] * _second_u[_qp].tr()) * this->tau() * U
    // *
    //             _grad_test[_i][_qp];
    if (_transient_term)
      residual += _rho[_qp] * _cp[_qp] * _u_dot[_qp] * this->tau() * U * _grad_test[_i][_qp];
  }

  return residual;
}

template class INSADTemperature<RESIDUAL>;
template class INSADTemperature<JACOBIAN>;
