//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#ifndef ADCOUPLEDCONVECTION_H_
#define ADCOUPLEDCONVECTION_H_

#include "ADKernel.h"

template <ComputeStage compute_stage>
class ADCoupledConvection;

template <>
InputParameters validParams<ADCoupledConvection<RESIDUAL>>();
template <>
InputParameters validParams<ADCoupledConvection<JACOBIAN>>();

/**
 * Define the ADKernel for a convection operator that looks like:
 *
 * grad_some_var dot u'
 *
 */
template <ComputeStage compute_stage>
class ADCoupledConvection : public ADKernel<compute_stage>
{
public:
  ADCoupledConvection(const InputParameters & parameters);

protected:
  virtual typename ResidualReturnType<compute_stage>::type computeQpResidual() override;

  using ADKernel<compute_stage>::_test;
  using ADKernel<compute_stage>::_i;
  using ADKernel<compute_stage>::_qp;
  using ADKernel<compute_stage>::_grad_u;

private:
  const typename VariableGradientType<compute_stage>::type & _velocity_vector;
};

#endif // ADCOUPLEDCONVECTION_H
