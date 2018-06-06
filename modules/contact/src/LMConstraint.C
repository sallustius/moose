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
#include "SystemBase.h"

registerMooseObject("ContactApp", LMConstraint);

template <>
InputParameters
validParams<LMConstraint>()
{
  InputParameters params = validParams<NodeFaceConstraint>();
  params.set<bool>("use_displaced_mesh") = true;

  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");

  return params;
}

LMConstraint::LMConstraint(const InputParameters & parameters)
  : NodeFaceConstraint(parameters), _disp_y_id(coupled("disp_y")), _disp_z_id(coupled("disp_z"))

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
  _connected_dof_indices.push_back(_var.nodalDofIndex());

  _qp = 0;

  _Kee(0, 0) += computeQpJacobian(Moose::SlaveSlave);
}

void
LMConstraint::computeOffDiagJacobian(unsigned jvar)
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
      Real a = -pinfo->_distance;
      Real b = _u_slave[_qp];
      return 1. - b / std::sqrt(a * a + b * b);
    }
  }
  return 0;
}

Real
LMConstraint::computeQpOffDiagJacobian(Moose::ConstraintJacobianType type, unsigned jvar)
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

      RealVectorValue distance_vec(*_current_node - pinfo->_closest_point);
      Real da_daj = 1. / -pinfo->_distance;

      if (jvar == _master_var_num)
        da_daj *= distance_vec(0);
      else if (jvar == _disp_y_id)
        da_daj *= distance_vec(1);
      else if (jvar == _disp_z_id)
        da_daj *= distance_vec(2);
      else
        da_daj *= 0;

      switch (type)
      {
        case Moose::SlaveSlave:
          da_daj *= 1;
          break;
        case Moose::SlaveMaster:
          da_daj *= -_phi_master[_j][_qp];
          break;
        default:
          mooseError("LMs do not have a master contribution.");
      }
      return da_daj - a / std::sqrt(a * a + b * b) * da_daj;
    }
  }
  return 0;
}
