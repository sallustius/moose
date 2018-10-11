//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADKERNEL_H
#define ADKERNEL_H

#include "KernelBase.h"

template <typename T, ComputeStage compute_stage>
class ADKernelTempl : public KernelBase, public MooseVariableInterface<T>
{
public:
  ADKernelTempl(const InputParameters & parameters);

  virtual ~ADKernelTempl();

  // See KernelBase base for documentation of these overridden methods
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  virtual MooseVariableFE<T> & variable() override { return _var; }

protected:
  /// Compute this Kernel's contribution to the residual at the current quadrature point
  virtual ADResidual computeQpResidual() = 0;

  /// This is a regular kernel so we cast to a regular MooseVariable
  MooseVariableFE<T> & _var;

  /// the current test function
  const ADTemplateVariableTestValue & _test;

  /// gradient of the test function
  const ADTemplateVariableTestGradient & _grad_test;

  /// Holds the solution at current quadrature points
  const ADTemplateVariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const ADTemplateVariableGradient & _grad_u;
};

template <ComputeStage compute_stage>
using ADKernel = ADKernelTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernel = ADKernelTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernel);
declareADValidParams(ADVectorKernel);

#define usingKernelMembers                                                                         \
  using ADKernel<compute_stage>::_test;                                                            \
  using ADKernel<compute_stage>::_qp;                                                              \
  using ADKernel<compute_stage>::_i;                                                               \
  using ADKernel<compute_stage>::_u;                                                               \
  using ADKernel<compute_stage>::_var;                                                             \
  using ADKernel<compute_stage>::_grad_test;                                                       \
  using ADKernel<compute_stage>::_grad_u

#define usingVectorKernelMembers                                                                   \
  using ADVectorKernel<compute_stage>::_test;                                                      \
  using ADVectorKernel<compute_stage>::_qp;                                                        \
  using ADVectorKernel<compute_stage>::_i;                                                         \
  using ADVectorKernel<compute_stage>::_u;                                                         \
  using ADVectorKernel<compute_stage>::_var;                                                       \
  using ADVectorKernel<compute_stage>::_grad_test;                                                 \
  using ADVectorKernel<compute_stage>::_grad_u

#endif /* ADKERNEL_H */
