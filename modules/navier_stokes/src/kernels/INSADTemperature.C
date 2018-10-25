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
                    ADKernel,
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
                    params.addParam<MaterialPropertyName>("cp_name", "cp", "specific heat name"););

template <ComputeStage compute_stage>
INSADTemperature<compute_stage>::INSADTemperature(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),

    // Coupled variables
    _u_vel(adCoupledValue("u")),
    _v_vel(adCoupledValue("v")),
    _w_vel(adCoupledValue("w")),

    // Material Properties
    _rho(adGetADMaterialProperty<Real>("rho_name")),
    _k(adGetADMaterialProperty<Real>("k_name")),
    _cp(adGetADMaterialProperty<Real>("cp_name"))
{
}

template <ComputeStage compute_stage>
ADResidual
INSADTemperature<compute_stage>::computeQpResidual()
{
  // The convection part, rho * cp u.grad(T) * v.
  // Note: _u is the temperature variable, _grad_u is its gradient.
  auto convective_part = _rho[_qp] * _cp[_qp] *
                         (_u_vel[_qp] * _grad_u[_qp](0) + _v_vel[_qp] * _grad_u[_qp](1) +
                          _w_vel[_qp] * _grad_u[_qp](2)) *
                         _test[_i][_qp];

  // Thermal conduction part, k * grad(T) * grad(v)
  auto conduction_part = _k[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];

  return convective_part + conduction_part;
}

template class INSADTemperature<RESIDUAL>;
template class INSADTemperature<JACOBIAN>;
