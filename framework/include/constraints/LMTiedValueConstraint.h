//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "NodeFaceConstraint.h"

// Forward Declarations
class LMTiedValueConstraint;

template <>
InputParameters validParams<LMTiedValueConstraint>();

/**
 * A LMTiedValueConstraint forces the value of a variable to be the same on both sides of an
 * interface.
 */
class LMTiedValueConstraint : public NodeFaceConstraint
{
public:
  static InputParameters validParams();

  LMTiedValueConstraint(const InputParameters & parameters);

  void computeResidual() override;
  void computeJacobian() override;
  void computeOffDiagJacobian(unsigned jvar) override;

protected:
  virtual Real computeQpSlaveValue() override;

  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::ConstraintJacobianType type, unsigned jvar) override;

  MooseVariable & _primal_var;
  const VariableValue & _primal_slave;
  const VariableValue & _primal_master;
  const unsigned int _primal_id;
};
