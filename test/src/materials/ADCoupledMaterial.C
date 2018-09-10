//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ADCoupledMaterial.h"

registerMooseObject("MooseTestApp", ADCoupledMaterial);

template <>
InputParameters
validParams<ADCoupledMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("coupled_var", "A coupledvariable");
  params.addRequiredParam<MaterialPropertyName>("ad_mat_prop",
                                                "Name of the ad property this material defines");
  params.addRequiredParam<MaterialPropertyName>(
      "regular_mat_prop", "Name of the regular property this material defines");
  return params;
}

ADCoupledMaterial::ADCoupledMaterial(const InputParameters & parameters)
  : Material(parameters),
    _ad_mat_prop(declareADProperty<Real>(getParam<MaterialPropertyName>("ad_mat_prop"))),
    _regular_mat_prop(declareProperty<Real>(getParam<MaterialPropertyName>("regular_mat_prop"))),
    _coupled_var(adCoupledValue("coupled_var"))
{
}

void
ADCoupledMaterial::computeQpProperties()
{
  _regular_mat_prop[_qp] = 4.0 * _coupled_var[_qp].value();
  _ad_mat_prop[_qp] = 4.0 * _coupled_var[_qp];
}
