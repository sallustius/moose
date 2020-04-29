//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADMatAdvection.h"

registerMooseObject("MooseApp", ADMatAdvection);

InputParameters
ADMatAdvection::validParams()
{
  auto params = ADKernelGrad::validParams();
  params.addClassDescription("Computes a conservative advection term (e.g. integrated by parts) in "
                             "which the advecting velocity is provided by a material");
  params.addRequiredParam<MaterialPropertyName>("vel", "advection velocity");
  params.addParam<MaterialPropertyName>(
      "advected_quantity",
      "An optional parameter to specify an advected quantity. If this parameter is not supplied, "
      "then it is assumed that the advected quantity is the variable that this kernel is acting "
      "on");
  return params;
}

ADMatAdvection::ADMatAdvection(const InputParameters & parameters)
  : ADKernelGrad(parameters),
    _vel(getADMaterialProperty<RealVectorValue>("vel")),
    _use_mat(isParamValid("advected_quantity")),
    _advected_quantity(_use_mat ? &getADMaterialProperty<Real>("advected_quantity") : nullptr)
{
}

ADRealVectorValue
ADMatAdvection::precomputeQpResidual()
{
  const ADReal & advected_quantity = _use_mat ? (*_advected_quantity)[_qp] : _u[_qp];

  return -_vel[_qp] * advected_quantity;
}
