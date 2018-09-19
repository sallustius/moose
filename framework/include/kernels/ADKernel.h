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

template <ComputeStage compute_stage>
class ADKernel;

template <>
InputParameters validParams<ADKernel<RESIDUAL>>();
template <>
InputParameters validParams<ADKernel<JACOBIAN>>();

template <ComputeStage compute_stage>
class ADKernel : public KernelBase, public MooseVariableInterface<Real>
{
public:
  ADKernel(const InputParameters & parameters);

  virtual ~ADKernel();

  // See KernelBase base for documentation of these overridden methods
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  virtual MooseVariable & variable() override { return _var; }

protected:
  /// Compute this Kernel's contribution to the residual at the current quadrature point
  virtual typename ResidualReturnType<compute_stage>::type computeQpResidual() = 0;

  /// This is a regular kernel so we cast to a regular MooseVariable
  MooseVariable & _var;

  /// the current test function
  const VariableTestValue & _test;

  /// gradient of the test function
  const VariableTestGradient & _grad_test;

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// Holds the solution at current quadrature points
  const typename VariableValueType<compute_stage>::type & _u;

  /// Holds the solution gradient at the current quadrature points
  const typename VariableGradientType<compute_stage>::type & _grad_u;

  /// Time derivative of u
  const VariableValue & _u_dot;

  /// Derivative of u_dot with respect to u
  const VariableValue & _du_dot_du;
};

#endif /* ADKERNEL_H */
