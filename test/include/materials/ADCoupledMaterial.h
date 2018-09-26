//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADCOUPLEDMATERIAL_H_
#define ADCOUPLEDMATERIAL_H_

#include "ADMaterial.h"

template <ComputeStage>
class ADCoupledMaterial;

template <>
InputParameters validParams<ADCoupledMaterial<RESIDUAL>>();
template <>
InputParameters validParams<ADCoupledMaterial<JACOBIAN>>();

/**
 * A material that couples a material property
 */
template <ComputeStage compute_stage>
class ADCoupledMaterial : public ADMaterial<compute_stage>
{
public:
  ADCoupledMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  typename MaterialPropertyType<compute_stage, Real>::type & _ad_mat_prop;
  MaterialProperty<Real> & _regular_mat_prop;

  const typename VariableValueType<compute_stage>::type & _coupled_var;

  using Material::_qp;
};

#endif // ADCOUPLEDMATERIAL_H
