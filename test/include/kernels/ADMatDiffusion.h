//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADMATDIFFUSION_H
#define ADMATDIFFUSION_H

#include "ADKernel.h"
#include "MaterialProperty.h"

// Forward Declaration
class ADMatDiffusion;

template <>
InputParameters validParams<ADMatDiffusion>();

class ADMatDiffusion : public ADKernel
{
public:
  ADMatDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual();

  const MaterialProperty<Real> & _regular_mat_prop;
  const MaterialProperty<ADReal> * _diff;
};

#endif // ADMATDIFFUSION_H
