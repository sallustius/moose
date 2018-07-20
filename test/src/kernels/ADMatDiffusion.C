//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ADMatDiffusion.h"

registerMooseObject("MooseTestApp", ADMatDiffusion);

template <>
InputParameters
validParams<ADMatDiffusion>()
{
  InputParameters params = validParams<ADKernel>();
  params.addRequiredParam<MaterialPropertyName>(
      "prop_name", "the name of the material property we are going to use");
  return params;
}

ADMatDiffusion::ADMatDiffusion(const InputParameters & parameters)
  : ADKernel(parameters),
    _regular_mat_prop(getMaterialProperty<Real>("mat_prop1")),
    _diff(getADMaterialProperty<Real>("prop_name"))
{
}

ADReal
ADMatDiffusion::computeQpResidual()
{
  return _diff[_qp] * _regular_mat_prop[_qp] * (_grad_test[_i][_qp] * _grad_u[_qp]);
}
