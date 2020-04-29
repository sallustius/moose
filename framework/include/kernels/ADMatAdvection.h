//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelGrad.h"

class ADMatAdvection : public ADKernelGrad
{
public:
  static InputParameters validParams();

  ADMatAdvection(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  /// The advecting velocity
  const ADMaterialProperty<RealVectorValue> & _vel;

  /// Whether to use a material property as the advected quantity
  const bool _use_mat;

  /// Advected material property. May not be used
  const ADMaterialProperty<Real> * _advected_quantity;
};
