//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NodalTangentialContactLM.h"
#include "PenetrationLocator.h"
#include "PenetrationInfo.h"
#include "SystemBase.h"
#include "Assembly.h"

#include "metaphysicl/dualnumberarray.h"
using MetaPhysicL::DualNumber;
using MetaPhysicL::NumberArray;

registerMooseObject("MooseTestApp", NodalTangentialContactLM);

template <>
InputParameters
validParams<NodalTangentialContactLM>()
{
  InputParameters params = validParams<NodeFaceConstraint>();
  params.set<bool>("use_displaced_mesh") = true;

  params.addRequiredCoupledVar(
      "contact_pressure",
      "The contact pressure. If using LM, this should be the normal lagrange multiplier");
  params.addRequiredCoupledVar("disp_y", "The y velocity");
  params.addRequiredParam<Real>("mu", "The coefficient of friction.");

  params.addClassDescription("Implements the KKT conditions for frictional Coulomb contact using "
                             "an NCP function. Requires that either the relative tangential "
                             "velocity is zero or the tangential stress is equal to the friction "
                             "coefficient times the normal contact pressure.");

  params.addRequiredParam<Real>("mu", "The friction coefficient for the Coulomb friction law");

  return params;
}

NodalTangentialContactLM::NodalTangentialContactLM(const InputParameters & parameters)
  : NodeFaceConstraint(parameters),
    _contact_pressure(getVar("contact_pressure", 0)->nodalValue()),
    _contact_pressure_id(coupled("contact_pressure")),
    _disp_x_dot(_master_var.nodalValueDot()),
    _disp_y_dot(getVar("disp_y", 0)->nodalValueDot()),
    _disp_y_id(coupled("disp_y")),
    _du_dot_du(_master_var.dofValuesDuDotDu()),

    _mu(getParam<Real>("mu")),
    _epsilon(std::numeric_limits<Real>::epsilon())
{
  _overwrite_slave_residual = false;
}

Real
NodalTangentialContactLM::computeQpSlaveValue()
{
  return _u_slave[_qp];
}

void
NodalTangentialContactLM::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());

  _qp = 0;

  re(0) = computeQpResidual(Moose::Slave);
}

void
NodalTangentialContactLM::computeJacobian()
{
  _Kee.resize(1, 1);
  // We have to calculate these connected dof indices because of logic in NonlinearSystemBase
  _connected_dof_indices.clear();
  _connected_dof_indices.push_back(_var.nodalDofIndex());

  _qp = 0;

  _Kee(0, 0) += computeQpJacobian(Moose::SlaveSlave);
}

void
NodalTangentialContactLM::computeOffDiagJacobian(unsigned jvar)
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
  _Kee(0, 0) = computeQpOffDiagJacobian(Moose::SlaveSlave, jvar);

  DenseMatrix<Number> & Ken =
      _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, _var.number(), jvar);

  auto master_jsize = var.dofIndicesNeighbor().size();

  for (_j = 0; _j < master_jsize; ++_j)
    Ken(0, _j) = computeQpOffDiagJacobian(Moose::SlaveMaster, jvar);
}

Real NodalTangentialContactLM::computeQpResidual(Moose::ConstraintType /*type*/)
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo)
    {
      RealVectorValue tangent_vec(pinfo->_normal(1), -pinfo->_normal(0), 0);
      auto v_dot_tan = RealVectorValue(_disp_x_dot, _disp_y_dot, 0) * tangent_vec;

      auto gap = -pinfo->_distance;

      return std::max(_mu * (_contact_pressure - gap), std::abs(_u_slave[_qp] + v_dot_tan)) *
                 _u_slave[_qp] -
             _mu * std::max(0., _contact_pressure - gap) * (_u_slave[_qp] + v_dot_tan);
    }
  }
  return _u_slave[_qp];
}

Real NodalTangentialContactLM::computeQpJacobian(Moose::ConstraintJacobianType /*type*/)
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo)
    {
      RealVectorValue tangent_vec(pinfo->_normal(1), -pinfo->_normal(0), 0);
      auto v_dot_tan = RealVectorValue(_disp_x_dot, _disp_y_dot, 0) * tangent_vec;

      auto gap = -pinfo->_distance;

      DualNumber<Real> dual_u_slave(_u_slave[_qp]);
      dual_u_slave.derivatives() = 1.;

      return (std::max(_mu * (_contact_pressure - gap), std::abs(dual_u_slave + v_dot_tan)) *
                  dual_u_slave -
              _mu * std::max(0., _contact_pressure - gap) * (dual_u_slave + v_dot_tan))
          .derivatives();
    }
  }
  return 1;
}

Real
NodalTangentialContactLM::computeQpOffDiagJacobian(Moose::ConstraintJacobianType type,
                                                   unsigned jvar)
{
  if (jvar == _master_var_num || jvar == _disp_y_id || jvar == _contact_pressure_id)
  {
    std::map<dof_id_type, PenetrationInfo *>::iterator found =
        _penetration_locator._penetration_info.find(_current_node->id());
    if (found != _penetration_locator._penetration_info.end())
    {
      PenetrationInfo * pinfo = found->second;
      if (pinfo)
      {
        // Our local dual number is going to depend on only three degrees of freedom: the slave
        // nodal dofs for disp_x (index 0), disp_y (index 1), and the contact pressure (index 2).
        // The latter of course exists only on the slave side
        typedef DualNumber<Real, NumberArray<3, Real>> LocalDN;

        RealVectorValue tangent_vec(pinfo->_normal(1), -pinfo->_normal(0), 0);

        LocalDN dual_disp_x_dot(_disp_x_dot);
        // We index the zeroth entry of the _du_dot_du member variable because there is only one
        // degree of freedom on the node
        dual_disp_x_dot.derivatives()[0] = _du_dot_du[0];

        LocalDN dual_disp_y_dot(_disp_y_dot);
        dual_disp_y_dot.derivatives()[1] = _du_dot_du[0];

        LocalDN dual_contact_pressure(_contact_pressure);
        dual_contact_pressure.derivatives()[2] = 1;

        auto v_dot_tan = VectorValue<LocalDN>(dual_disp_x_dot, dual_disp_y_dot, 0) * tangent_vec;

        LocalDN gap(-pinfo->_distance);

        if (jvar == _master_var_num)
        {
          switch (type)
          {
            case Moose::SlaveSlave:
              gap.derivatives()[0] = pinfo->_normal(0);
              break;
            case Moose::SlaveMaster:
              gap.derivatives()[0] = pinfo->_normal(0) * -_phi_master[_j][_qp];
              break;
            default:
              mooseError(
                  "You shouldn't be calling in with types other than SlaveSlave or SlaveMaster");
              break;
          }
        }
        else if (jvar == _disp_y_id)
        {
          switch (type)
          {
            case Moose::SlaveSlave:
              gap.derivatives()[1] = pinfo->_normal(1);
              break;
            case Moose::SlaveMaster:
              gap.derivatives()[1] = pinfo->_normal(1) * -_phi_master[_j][_qp];
              break;
            default:
              mooseError(
                  "You shouldn't be calling in with types other than SlaveSlave or SlaveMaster");
              break;
          }
        }

        auto ncp_value =
            std::max(_mu * (dual_contact_pressure - gap), std::abs(_u_slave[_qp] + v_dot_tan)) *
                _u_slave[_qp] -
            _mu * std::max(0., dual_contact_pressure - gap) * (_u_slave[_qp] + v_dot_tan);

        if (jvar == _contact_pressure_id)
          return ncp_value.derivatives()[2];
        else if (jvar == _master_var.number())
          return ncp_value.derivatives()[0];
        else if (jvar == _disp_y_id)
          return ncp_value.derivatives()[1];
      }
    }
  }
  return 0.;
}
