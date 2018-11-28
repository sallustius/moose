//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADBase.h"
#include "Assembly.h"

defineADValidParams(
    INSADBase,
    ADKernel,
    params.addClassDescription("This class computes various strong and weak components of the "
                               "incompressible navier stokes equations which can then be assembled "
                               "together in child classes.");
    // Coupled variables
    params.addRequiredCoupledVar("u", "x-velocity");
    params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
    params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
    params.addRequiredCoupledVar("p", "pressure");

    params.addParam<RealVectorValue>("gravity",
                                     RealVectorValue(0, 0, 0),
                                     "Direction of the gravity vector");

    params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
    params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");

    params.addParam<Real>("alpha", 1., "Multiplicative factor on the stabilization parameter tau.");
    params.addParam<bool>("laplace",
                          true,
                          "Whether the viscous term of the momentum equations is in laplace form.");
    params.addParam<bool>("convective_term", true, "Whether to include the convective term.");
    params.addParam<bool>("transient_term",
                          false,
                          "Whether there should be a transient term in the momentum residuals."););

template <ComputeStage compute_stage>
INSADBase<compute_stage>::INSADBase(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),

    // Coupled variables
    _u_vel(adCoupledValue("u")),
    _v_vel(adCoupledValue("v")),
    _w_vel(adCoupledValue("w")),
    _p(adCoupledValue("p")),

    // Gradients
    _grad_u_vel(adCoupledGradient("u")),
    _grad_v_vel(adCoupledGradient("v")),
    _grad_w_vel(adCoupledGradient("w")),
    _grad_p(adCoupledGradient("p")),

    // second derivative tensors
    _second_u_vel(adCoupledSecond("u")),
    _second_v_vel(adCoupledSecond("v")),
    _second_w_vel(adCoupledSecond("w")),

    // time derivatives
    _u_vel_dot(adCoupledDot("u")),
    _v_vel_dot(adCoupledDot("v")),
    _w_vel_dot(adCoupledDot("w")),

    _gravity(adGetParam<RealVectorValue>("gravity")),

    // Material properties
    _mu(adGetADMaterialProperty<Real>("mu_name")),
    _rho(adGetADMaterialProperty<Real>("rho_name")),

    _alpha(adGetParam<Real>("alpha")),
    _laplace(adGetParam<bool>("laplace")),
    _convective_term(adGetParam<bool>("convective_term")),
    _transient_term(adGetParam<bool>("transient_term"))
{
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::convectiveTerm()
{
  INSVectorValue<compute_stage> U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  return _rho[_qp] * INSVectorValue<compute_stage>(
                         U * _grad_u_vel[_qp], U * _grad_v_vel[_qp], U * _grad_w_vel[_qp]);
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::strongViscousTermLaplace()
{
  return -_mu[_qp] * INSVectorValue<compute_stage>(
                         _second_u_vel[_qp].tr(), _second_v_vel[_qp].tr(), _second_w_vel[_qp].tr());
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::strongViscousTermTraction()
{
  return strongViscousTermLaplace() -
         _mu[_qp] *
             (_second_u_vel[_qp].row(0) + _second_v_vel[_qp].row(1) + _second_w_vel[_qp].row(2));
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::weakViscousTermLaplace(unsigned comp)
{
  switch (comp)
  {
    case 0:
      return _mu[_qp] * _grad_u_vel[_qp];

    case 1:
      return _mu[_qp] * _grad_v_vel[_qp];

    case 2:
      return _mu[_qp] * _grad_w_vel[_qp];

    default:
      return adZeroValue()[_qp];
  }
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::weakViscousTermTraction(unsigned comp)
{
  switch (comp)
  {
    case 0:
    {
      INSVectorValue<compute_stage> transpose(
          _grad_u_vel[_qp](0), _grad_v_vel[_qp](0), _grad_w_vel[_qp](0));
      return _mu[_qp] * _grad_u_vel[_qp] + _mu[_qp] * transpose;
    }

    case 1:
    {
      INSVectorValue<compute_stage> transpose(
          _grad_u_vel[_qp](1), _grad_v_vel[_qp](1), _grad_w_vel[_qp](1));
      return _mu[_qp] * _grad_v_vel[_qp] + _mu[_qp] * transpose;
    }

    case 2:
    {
      INSVectorValue<compute_stage> transpose(
          _grad_u_vel[_qp](2), _grad_v_vel[_qp](2), _grad_w_vel[_qp](2));
      return _mu[_qp] * _grad_w_vel[_qp] + _mu[_qp] * transpose;
    }

    default:
      return adZeroValue()[_qp];
  }
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::strongPressureTerm()
{
  return _grad_p[_qp];
}

template <ComputeStage compute_stage>
INSReal<compute_stage>
INSADBase<compute_stage>::weakPressureTerm()
{
  return -_p[_qp];
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::gravityTerm()
{
  return -_rho[_qp] * _gravity;
}

template <ComputeStage compute_stage>
INSVectorValue<compute_stage>
INSADBase<compute_stage>::timeDerivativeTerm()
{
  return _rho[_qp] *
         INSVectorValue<compute_stage>(_u_vel_dot[_qp], _v_vel_dot[_qp], _w_vel_dot[_qp]);
}

template <ComputeStage compute_stage>
INSReal<compute_stage>
INSADBase<compute_stage>::tau()
{
  auto && nu = _mu[_qp] / _rho[_qp];
  INSVectorValue<compute_stage> U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  auto && h = _current_elem->hmax();
  auto && transient_part = _transient_term ? 4. / (_dt * _dt) : 0.;
  return _alpha / std::sqrt(transient_part + (2. * U.norm() / h) * (2. * U.norm() / h) +
                            9. * (4. * nu / (h * h)) * (4. * nu / (h * h)));
}

template <ComputeStage compute_stage>
INSReal<compute_stage>
INSADBase<compute_stage>::tauNodal()
{
  auto nu = _mu[_qp] / _rho[_qp];
  INSVectorValue<compute_stage> U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  const auto & h = _current_elem->hmax();
  INSReal<compute_stage> xi;
  if (nu < std::numeric_limits<Real>::epsilon())
    xi = 1;
  else
  {
    auto alpha = U.norm() * h / (2 * nu);
    xi = 1. / std::tanh(alpha) - 1. / alpha;
  }
  return h / (2 * U.norm()) * xi;
}

template class INSADBase<RESIDUAL>;
template class INSADBase<JACOBIAN>;
