//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADMATTYPES_H_
#define ADMATTYPES_H_

#include "Material.h"

class ADMatTypes;

template <>
InputParameters validParams<ADMatTypes>();

/**
 * A material that couples a material property
 */
class ADMatTypes : public Material
{
public:
  ADMatTypes(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  ADMaterialProperty<Real> & _scalar_ad_prop;
  ADMaterialProperty<RealVectorValue> & _vector_ad_prop;
  ADMaterialProperty<RealTensorValue> & _tensor_ad_prop;
  MaterialProperty<Real> & _scalar_reg_prop;
  MaterialProperty<RealVectorValue> & _vector_reg_prop;
  MaterialProperty<RealTensorValue> & _tensor_reg_prop;
};

#endif // ADMATTYPES_H
