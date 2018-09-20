//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADVALUETEST_H_
#define ADVALUETEST_H_

#include "ADKernel.h"

template <ComputeStage>
class ADValueTest;

template <>
InputParameters validParams<ADValueTest<RESIDUAL>>();
template <>
InputParameters validParams<ADValueTest<JACOBIAN>>();

template <ComputeStage compute_stage>
class ADValueTest : public ADKernel<compute_stage>
{
public:
  ADValueTest(const InputParameters & parameters);

protected:
  virtual typename ResidualReturnType<compute_stage>::type computeQpResidual();

  using ADKernel<compute_stage>::_u;
  using ADKernel<compute_stage>::_test;
  using ADKernel<compute_stage>::_i;
  using ADKernel<compute_stage>::_qp;
};

#endif /* ADVALUETEST_H_ */
