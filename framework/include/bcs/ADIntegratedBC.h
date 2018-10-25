//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADINTEGRATEDBC_H
#define ADINTEGRATEDBC_H

#include "IntegratedBCBase.h"
#include "MooseVariableInterface.h"

/**
 * Base class for deriving any boundary condition of a integrated type
 */
template <typename T, ComputeStage compute_stage>
class ADIntegratedBCTempl : public IntegratedBCBase, public MooseVariableInterface<T>
{
public:
  ADIntegratedBCTempl(const InputParameters & parameters);

  virtual MooseVariableFE<T> & variable() override { return _var; }

  void computeResidual() override;
  void computeJacobian() override;
  void computeJacobianBlock(MooseVariableFEBase & jvar) override;
  void computeJacobianBlockScalar(unsigned int jvar) override;

protected:
  /**
   * Compute this IntegratedBC's contribution to the residual at the current quadrature point
   */
  virtual ADResidual computeQpResidual() = 0;

  /// The variable that this IntegratedBC operates on
  MooseVariableFE<T> & _var;

  /// normals at quadrature points
  const MooseArray<Point> & _normals;

  // test functions

  /// test function values (in QPs)
  const ADTemplateVariableTestValue & _test;
  /// gradients of test functions  (in QPs)
  const ADTemplateVariableTestGradient & _grad_test;

  /// the values of the unknown variable this BC is acting on
  const ADTemplateVariableValue & _u;
  /// the gradient of the unknown variable this BC is acting on
  const ADTemplateVariableGradient & _grad_u;
};

template <ComputeStage compute_stage>
using ADIntegratedBC = ADIntegratedBCTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorIntegratedBC = ADIntegratedBCTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADIntegratedBC);
declareADValidParams(ADVectorIntegratedBC);

#define usingIntegratedBCMembers                                                                   \
  using ADIntegratedBC<compute_stage>::_test;                                                      \
  using ADIntegratedBC<compute_stage>::_qp;                                                        \
  using ADIntegratedBC<compute_stage>::_i;                                                         \
  using ADIntegratedBC<compute_stage>::_u;                                                         \
  using ADIntegratedBC<compute_stage>::_var;                                                       \
  using ADIntegratedBC<compute_stage>::_grad_test;                                                 \
  using ADIntegratedBC<compute_stage>::_grad_u

#define usingVectorIntegratedBCMembers usingIntegratedBCMembers

#endif /* ADINTEGRATEDBC_H */
