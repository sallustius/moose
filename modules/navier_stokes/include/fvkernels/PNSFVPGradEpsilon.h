//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVElementalKernel.h"

class Function;

class PNSFVPGradEpsilon : public FVElementalKernel
{
public:
  static InputParameters validParams();
  PNSFVPGradEpsilon(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// The gradient of the pressure
  const ADMaterialProperty<Real> & _pressure;
  const Function & _eps_function;

  /// index x|y|z
  const unsigned int _index;
};
