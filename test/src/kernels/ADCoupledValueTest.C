//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ADCoupledValueTest.h"

registerADMooseObject("MooseTestApp", ADCoupledValueTest);

template <>
InputParameters
validParams<ADCoupledValueTest<RESIDUAL>>()
{
  InputParameters params = validParams<ADKernel<RESIDUAL>>();
  params.addCoupledVar("v", 2.0, "The coupled variable.");
  return params;
}
template <>
InputParameters
validParams<ADCoupledValueTest<JACOBIAN>>()
{
  return validParams<ADCoupledValueTest<RESIDUAL>>();
}

template <ComputeStage compute_stage>
ADCoupledValueTest<compute_stage>::ADCoupledValueTest(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters), _v(this->template adCoupledValue<compute_stage>("v"))
{
}

template <ComputeStage compute_stage>
typename ResidualReturnType<compute_stage>::type
ADCoupledValueTest<compute_stage>::computeQpResidual()
{
  return _test[_i][_qp] * -_v[_qp];
}
