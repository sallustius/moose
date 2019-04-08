//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BadTimeDerivative.h"
#include "Function.h"

registerADMooseObject("NavierStokesTestApp", BadTimeDerivative);

defineADValidParams(BadTimeDerivative, ADKernel, );

template <ComputeStage compute_stage>
BadTimeDerivative<compute_stage>::BadTimeDerivative(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters), _u_dot(_var.template adUDot<compute_stage>())
{
}

template <ComputeStage compute_stage>
ADReal
BadTimeDerivative<compute_stage>::computeQpResidual()
{
  return _test[_i][_qp] * _u_dot[_qp];
}
