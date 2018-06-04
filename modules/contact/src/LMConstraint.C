//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LMConstraint.h"
#include "PenetrationLocator.h"
#include "PenetrationInfo.h"

registerMooseObject("ContactApp", LMConstraint);

template <>
InputParameters
validParams<LMConstraint>()
{
  InputParameters params = validParams<NodeFaceConstraint>();
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

LMConstraint::LMConstraint(const InputParameters & parameters) : NodeFaceConstraint(parameters)
{
  _overwrite_slave_residual = false;
}

Real
LMConstraint::computeQpSlaveValue()
{
  return _u_slave[_qp];
}

void
LMConstraint::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());

  _qp = 0;

  re(0) = computeQpResidual(Moose::Slave);
}

void
LMConstraint::computeJacobian()
{
  _Kee.resize(1, 1);
  _connected_dof_indices.clear();
  _connected_dof_indices.push_back(variable().nodalDofIndex());

  _qp = 0;

  _Kee(0, 0) += computeQpJacobian(Moose::SlaveSlave);
}

Real LMConstraint::computeQpResidual(Moose::ConstraintType /*type*/)
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo != NULL)
    {
      Real a = -pinfo->_distance;
      Real b = _u_slave[_qp];
      return a + b - std::sqrt(a * a + b * b);
    }
  }
  return 0;
}

Real LMConstraint::computeQpJacobian(Moose::ConstraintJacobianType /*type*/)
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo != NULL)
    {
      Real a = pinfo->_distance;
      Real b = _u_slave[_qp];
      return 1. + b / std::sqrt(a * a + b * b);
    }
  }
  return 0;
}
