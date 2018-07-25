//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ADMatTypes.h"

registerMooseObject("MooseTestApp", ADMatTypes);

template <>
InputParameters
validParams<ADMatTypes>()
{
  return validParams<Material>();
}

ADMatTypes::ADMatTypes(const InputParameters & parameters)
  : Material(parameters),
    _scalar_ad_prop(declareADProperty<Real>("scalar_ad_prop")),
    _vector_ad_prop(declareADProperty<RealVectorValue>("vector_ad_prop")),
    _tensor_ad_prop(declareADProperty<RealTensorValue>("tensor_ad_prop")),
    _scalar_reg_prop(declareProperty<Real>("scalar_reg_prop")),
    _vector_reg_prop(declareProperty<RealVectorValue>("vector_reg_prop")),
    _tensor_reg_prop(declareProperty<RealTensorValue>("tensor_reg_prop"))
{
}

void
ADMatTypes::computeQpProperties()
{
  // All binary and unary operators are implemented using the same macro, so we're just going to use
  // multiply (*) as a surrogate to represent them all
  _scalar_ad_prop[_qp] = 1.;
  _scalar_ad_prop[_qp] *= _scalar_ad_prop[_qp];
  _scalar_ad_prop[_qp] = _scalar_ad_prop[_qp] * _scalar_ad_prop[_qp];

  _vector_ad_prop[_qp] = _scalar_ad_prop[_qp] * RealVectorValue(1., 1., 1.);
  _scalar_ad_prop[_qp] = operator*<Real>(_vector_ad_prop[_qp], _vector_ad_prop[_qp]);
  // _vector_ad_prop[_qp] *= 2.;
  // _vector_ad_prop[_qp] = 2. * _vector_ad_prop[_qp] + _vector_ad_prop[_qp] * 2.;

  // _tensor_ad_prop[_qp] = _scalar_ad_prop[_qp] * RealTensorValue(1., 1., 1., 1., 1., 1., 1.,
  // 1., 1.); _vector_ad_prop[_qp] = _tensor_ad_prop[_qp] * _vector_ad_prop[_qp];
  // _tensor_ad_prop[_qp] *= 2.;
  // _tensor_ad_prop[_qp] *= _tensor_ad_prop[_qp];
  // _tensor_ad_prop[_qp] = _tensor_ad_prop[_qp] * _tensor_ad_prop[_qp];
  // _tensor_ad_prop[_qp] = 2. * _tensor_ad_prop[_qp] + _tensor_ad_prop[_qp] * 2.;

  // _scalar_reg_prop[_qp] = _scalar_ad_prop[_qp];
  // _vector_reg_prop[_qp] = _vector_ad_prop[_qp];
  // _tensor_reg_prop[_qp] = _tensor_ad_prop[_qp];
}
