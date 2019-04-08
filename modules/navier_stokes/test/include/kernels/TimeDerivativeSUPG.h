//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef TIMEDERIVATIVESUPG_H
#define TIMEDERIVATIVESUPG_H

#include "ADKernel.h"

// Forward Declarations
template <ComputeStage>
class TimeDerivativeSUPG;

declareADValidParams(TimeDerivativeSUPG);

template <ComputeStage compute_stage>
class TimeDerivativeSUPG : public ADKernel<compute_stage>
{
public:
  TimeDerivativeSUPG(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const RealVectorValue _velocity;
  const ADVariableValue & _u_dot;
  const Real _diff;

  usingKernelMembers;
};

#endif
