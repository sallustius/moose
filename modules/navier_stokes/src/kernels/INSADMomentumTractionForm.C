//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADMomentumTractionForm.h"

registerADMooseObject("NavierStokesApp", INSADMomentumTractionForm);

defineADValidParams(INSADMomentumTractionForm,
                    INSADMomentumBase,
                    params.addClassDescription(
                        "This class computes momentum equation residual and Jacobian viscous "
                        "contributions for the 'traction' form of the governing equations."););

template <ComputeStage compute_stage>
INSADMomentumTractionForm<compute_stage>::INSADMomentumTractionForm(
    const InputParameters & parameters)
  : INSADMomentumBase<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
ADResidual
INSADMomentumTractionForm<compute_stage>::computeQpResidualViscousPart()
{
  // The component'th row (or col, it's symmetric) of the viscous stress tensor
  INSVectorValue<compute_stage> tau_row;

  switch (_component)
  {
    case 0:
      tau_row(0) = 2. * _grad_u_vel[_qp](0);                  // 2*du/dx1
      tau_row(1) = _grad_u_vel[_qp](1) + _grad_v_vel[_qp](0); // du/dx2 + dv/dx1
      tau_row(2) = _grad_u_vel[_qp](2) + _grad_w_vel[_qp](0); // du/dx3 + dw/dx1
      break;

    case 1:
      tau_row(0) = _grad_v_vel[_qp](0) + _grad_u_vel[_qp](1); // dv/dx1 + du/dx2
      tau_row(1) = 2. * _grad_v_vel[_qp](1);                  // 2*dv/dx2
      tau_row(2) = _grad_v_vel[_qp](2) + _grad_w_vel[_qp](1); // dv/dx3 + dw/dx2
      break;

    case 2:
      tau_row(0) = _grad_w_vel[_qp](0) + _grad_u_vel[_qp](2); // dw/dx1 + du/dx3
      tau_row(1) = _grad_w_vel[_qp](1) + _grad_v_vel[_qp](2); // dw/dx2 + dv/dx3
      tau_row(2) = 2. * _grad_w_vel[_qp](2);                  // 2*dw/dx3
      break;

    default:
      mooseError("Unrecognized _component requested.");
  }

  // The viscous part, _mu[_qp] * tau : grad(v)
  return _mu[_qp] * (tau_row * _grad_test[_i][_qp]);
}

template class INSADMomentumTractionForm<RESIDUAL>;
template class INSADMomentumTractionForm<JACOBIAN>;
