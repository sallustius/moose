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
  virtual void computeADOffDiagJacobian() override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;
  virtual void beforeTestLoop() {}
  virtual void beforeQpLoop() {}

  virtual MooseVariableFE<T> & variable() override { return _var; }

protected:
  /// Compute this Kernel's contribution to the residual at the current quadrature point
  virtual ADResidual computeQpResidual() = 0;

  /// This is a regular kernel so we cast to a regular MooseVariable
  MooseVariableFE<T> & _var;

  /// the current test function
  const ADTemplateVariableTestValue & _test;

  /// gradient of the test function
  const typename VariableTestGradientType<compute_stage, T>::type & _grad_test;

  /// Holds the solution at current quadrature points
  const ADTemplateVariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const ADTemplateVariableGradient & _grad_u;

  /// The ad version of JxW
  const MooseArray<typename Moose::RealType<compute_stage>::type> & _ad_JxW;
};

template <ComputeStage compute_stage>
using ADKernel = ADKernelTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernel = ADKernelTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernel);
declareADValidParams(ADVectorKernel);

#define usingTemplKernelMembers(type)                                                              \
  using ADKernelTempl<type, compute_stage>::_test;                                                 \
  using ADKernelTempl<type, compute_stage>::_qp;                                                   \
  using ADKernelTempl<type, compute_stage>::_i;                                                    \
  using ADKernelTempl<type, compute_stage>::_u;                                                    \
  using ADKernelTempl<type, compute_stage>::_var;                                                  \
  using ADKernelTempl<type, compute_stage>::_grad_test;                                            \
  using ADKernelTempl<type, compute_stage>::_grad_u;                                               \
  using ADKernelTempl<type, compute_stage>::_dt;                                                   \
  using ADKernelTempl<type, compute_stage>::_current_elem;                                         \
  using ADKernelTempl<type, compute_stage>::_t;                                                    \
  using ADKernelTempl<type, compute_stage>::_q_point;                                              \
  using ADKernelTempl<type, compute_stage>::_displacements;                                        \
  using ADKernelTempl<type, compute_stage>::beforeTestLoop;                                        \
  using ADKernelTempl<type, compute_stage>::beforeQpLoop;                                          \
  using ADKernelTempl<type, compute_stage>::getFunction

#define usingKernelMembers usingTemplKernelMembers(Real)
#define usingVectorKernelMembers usingTemplKernelMembers(RealVectorValue)

#endif /* ADKERNEL_H */
