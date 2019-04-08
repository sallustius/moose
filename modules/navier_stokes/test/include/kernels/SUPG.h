//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SUPG_H
#define SUPG_H

#include "ADKernel.h"

// Forward Declarations
template <ComputeStage>
class SUPG;

declareADValidParams(SUPG);

template <ComputeStage compute_stage>
class SUPG : public ADKernel<compute_stage>
{
public:
  SUPG(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  Function & _ffn;
  const RealVectorValue _velocity;
  const ADVariableSecond & _second_u;
  const ADVariableValue & _u_dot;
  const Real _diff;
  const bool _include_transient_term;

  MooseEnum _tau_type;

  using KernelBase::_q_point;

  usingKernelMembers;
};

#endif
