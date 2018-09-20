//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADCOUPLEDVALUETEST_H_
#define ADCOUPLEDVALUETEST_H_

#include "ADKernel.h"

template <ComputeStage compute_stage>
class ADCoupledValueTest : public ADKernel<compute_stage>
{
public:
  ADCoupledValueTest(const InputParameters & parameters);

protected:
  virtual typename ResidualReturnType<compute_stage>::type computeQpResidual();

  const typename VariableValueType<compute_stage>::type & _v;

  using ADKernel<compute_stage>::_test;
  using ADKernel<compute_stage>::_i;
  using ADKernel<compute_stage>::_qp;
  using Coupleable::adCoupledValue;
};

template <>
InputParameters validParams<ADCoupledValueTest<RESIDUAL>>();
template <>
InputParameters validParams<ADCoupledValueTest<JACOBIAN>>();

#endif /* ADCOUPLEDVALUETEST_H_ */
