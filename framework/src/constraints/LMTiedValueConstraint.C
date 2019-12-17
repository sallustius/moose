//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LMTiedValueConstraint.h"

// MOOSE includes
#include "MooseVariableFE.h"
#include "SystemBase.h"

#include "libmesh/sparse_matrix.h"

registerMooseObject("MooseApp", LMTiedValueConstraint);

defineLegacyParams(LMTiedValueConstraint);

InputParameters
LMTiedValueConstraint::validParams()
{
  InputParameters params = NodeFaceConstraint::validParams();
  params.set<bool>("use_displaced_mesh") = true;
  params.addRequiredCoupledVar("primal_var", "the primal variable");
  return params;
}

LMTiedValueConstraint::LMTiedValueConstraint(const InputParameters & parameters)
  : NodeFaceConstraint(parameters),
    _primal_var(*getVar("primal_var", 0)),
    _primal_slave(_primal_var.dofValues()),
    _primal_master(_primal_var.slnNeighbor()),
    _primal_id(_primal_var.number())
{
  _overwrite_slave_residual = false;
}

Real
LMTiedValueConstraint::computeQpSlaveValue()
{
  return _primal_master[_qp];
}

void
LMTiedValueConstraint::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());

  _qp = 0;

  _slave_residual = re(0) = computeQpResidual(Moose::Slave);

  _slave_residual_computed = true;
}

void
LMTiedValueConstraint::computeJacobian()
{
  _Kee.resize(1, 1);
  _connected_dof_indices.clear();
  _connected_dof_indices.push_back(_var.nodalDofIndex());

  _qp = 0;

  _Kee(0, 0) += computeQpJacobian(Moose::SlaveSlave);
}

void
LMTiedValueConstraint::computeOffDiagJacobian(unsigned jvar)
{
  if (jvar == _var.number())
  {
    computeJacobian();
    return;
  }

  MooseVariableFEBase & var = _sys.getVariable(0, jvar);
  _connected_dof_indices.clear();
  _connected_dof_indices.push_back(var.nodalDofIndex());

  _qp = 0;

  _Kee.resize(1, 1);
  _Kee(0, 0) += computeQpOffDiagJacobian(Moose::SlaveSlave, jvar);

  DenseMatrix<Number> & Ken =
      _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, _var.number(), jvar);

  for (_j = 0; _j < _phi_master.size(); ++_j)
    Ken(0, _j) += computeQpOffDiagJacobian(Moose::SlaveMaster, jvar);
}

Real
LMTiedValueConstraint::computeQpResidual(Moose::ConstraintType type)
{
  Real retVal = 0;
  switch (type)
  {
    case Moose::Slave:
      retVal = _primal_slave[_qp] - _primal_master[_qp];
      break;
    default:
      break;
  }
  return retVal;
}

Real LMTiedValueConstraint::computeQpJacobian(Moose::ConstraintJacobianType) { return 0; }

Real
LMTiedValueConstraint::computeQpOffDiagJacobian(Moose::ConstraintJacobianType type, unsigned jvar)
{
  if (jvar != _primal_var.number())
    return 0;

  switch (type)
  {
    case Moose::SlaveSlave:
      return 1;
    case Moose::SlaveMaster:
      return -_phi_master[_j][_qp];
    default:
      return 0;
  }
}
