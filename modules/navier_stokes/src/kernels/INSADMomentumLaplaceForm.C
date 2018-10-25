//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADMomentumLaplaceForm.h"

registerADMooseObject("NavierStokesApp", INSADMomentumLaplaceForm);

defineADValidParams(INSADMomentumLaplaceForm,
                    INSADMomentumBase,
                    params.addClassDescription(
                        "This class computes momentum equation residual and Jacobian viscous "
                        "contributions for the 'Laplacian' form of the governing equations."););

template <ComputeStage compute_stage>
INSADMomentumLaplaceForm<compute_stage>::INSADMomentumLaplaceForm(
    const InputParameters & parameters)
  : INSADMomentumBase<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
ADResidual
INSADMomentumLaplaceForm<compute_stage>::computeQpResidualViscousPart()
{
  // Simplified version: mu * Laplacian(u_component)
  return _mu[_qp] * (_grad_u[_qp] * _grad_test[_i][_qp]);
}

template class INSADMomentumLaplaceForm<RESIDUAL>;
template class INSADMomentumLaplaceForm<JACOBIAN>;
