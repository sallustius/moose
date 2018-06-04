//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef LMCONSTRAINT_H
#define LMCONSTRAINT_H

// MOOSE includes
#include "NodeFaceConstraint.h"

// Forward Declarations
class LMConstraint;

template <>
InputParameters validParams<LMConstraint>();

class LMConstraint : public NodeFaceConstraint
{
public:
  LMConstraint(const InputParameters & parameters);

protected:
  virtual Real computeQpSlaveValue() override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;
};

#endif
