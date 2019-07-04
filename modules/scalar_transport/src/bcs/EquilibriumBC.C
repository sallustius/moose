//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EquilibriumBC.h"

registerADMooseObject("ScalarTransportApp", EquilibriumBC);

defineADValidParams(EquilibriumBC, ADIntegratedBC, );

template <ComputeStage compute_stage>
EquilibriumBC<compute_stage>::EquilibriumBC(const InputParameters & parameters)
  : ADIntegratedBC<compute_stage>(parameters)
{
}

template <ComputeStage compute_stage>
ADReal
EquilibriumBC<compute_stage>::computeQpResidual()
{
  return _test[_i][_qp] * 2. * _u[_qp];
}
