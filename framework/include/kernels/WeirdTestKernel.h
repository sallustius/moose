//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef WEIRDTESTKERNEL_H
#define WEIRDTESTKERNEL_H

#include "ADKernel.h"

template <ComputeStage>
class WeirdTestKernel;

declareADValidParams(WeirdTestKernel);

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
template <ComputeStage compute_stage>
class WeirdTestKernel : public ADKernel<compute_stage>
{
public:
  WeirdTestKernel(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  usingKernelMembers;
};

#endif /* WEIRDTESTKERNEL_H */
