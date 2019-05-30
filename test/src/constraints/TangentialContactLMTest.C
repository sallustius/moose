//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TangentialContactLMTest.h"

registerADMooseObject("MooseTestApp", TangentialContactLMTest);

defineADValidParams(
    TangentialContactLMTest,
    ADMortarConstraint,
    params.addRequiredParam<NonlinearVariableName>("slave_disp_y",
                                                   "The y displacement variable on the slave face");
    params.addParam<NonlinearVariableName>("master_disp_y",
                                           "The y displacement variable on the master face");
    params.addRequiredParam<NonlinearVariableName>(
        "contact_pressure",
        "The normal contact pressure; oftentimes this may be a separate lagrange multiplier "
        "variable");
    params.addRequiredParam<Real>("friction_coefficient", "The friction coefficient");

    MooseEnum ncp_type("min fb", "min");
    params.addParam<MooseEnum>("ncp_function_type",
                               ncp_type,
                               "The type of the nonlinear complimentarity function; options are "
                               "min or fb where fb stands for Fischer-Burmeister"););

template <ComputeStage compute_stage>
TangentialContactLMTest<compute_stage>::TangentialContactLMTest(const InputParameters & parameters)
  : ADMortarConstraint<compute_stage>(parameters),
    _slave_disp_y(_subproblem.getStandardVariable(_tid, parameters.getMooseType("slave_disp_y"))),
    _master_disp_y(
        isParamValid("master_disp_y")
            ? _subproblem.getStandardVariable(_tid, parameters.getMooseType("master_disp_y"))
            : _subproblem.getStandardVariable(_tid, parameters.getMooseType("slave_disp_y"))),
    _contact_pressure_var(
        _subproblem.getStandardVariable(_tid, parameters.getMooseType("contact_pressure"))),
    _contact_pressure(_contact_pressure_var.template adSlnLower<compute_stage>()),
    _slave_x_dot(_slave_var.template adUDot<compute_stage>()),
    _master_x_dot(_master_var.template adUDotNeighbor<compute_stage>()),
    _slave_y_dot(_slave_disp_y.template adUDot<compute_stage>()),
    _master_y_dot(_master_disp_y.template adUDotNeighbor<compute_stage>()),
    _friction_coeff(adGetParam<Real>("friction_coefficient")),
    _epsilon(std::numeric_limits<Real>::epsilon()),
    _ncp_type(adGetParam<MooseEnum>("ncp_function_type"))
{
}

template <ComputeStage compute_stage>
ADReal
TangentialContactLMTest<compute_stage>::computeQpResidual(Moose::MortarType mortar_type)
{
  switch (mortar_type)
  {
    case Moose::MortarType::Lower:
    {
      // Check whether we project onto a master face
      if (_has_master)
      {
        // Build the velocity vector
        ADRealVectorValue relative_velocity(
            _slave_x_dot[_qp] - _master_x_dot[_qp], _slave_y_dot[_qp] - _master_y_dot[_qp], 0);

        // Get the component in the tangential direction
        auto v_dot_tan = relative_velocity * _tangents[_qp][0];

        auto gap_vec = gapVec<compute_stage>();
        auto gap = gap_vec * _normals[_qp];

        return _test[_i][_qp] *
                   std::max(_mu * (_contact_pressure[_qp] - gap),
                            std::abs(_lambda[_qp] + v_dot_tan)) *
                   _lambda[_qp] -
               _mu * std::max(0., _contact_pressure[_qp] - gap) * (_lambda[_qp] + v_dot_tan);
      }
      else
        // If not in contact then we force the tangential lagrange multiplier to zero
        return _test[_i][_qp] * _lambda[_qp];
    }

    default:
      return 0;
  }
}

template <>
RealVectorValue
TangentialContactLMTest<RESIDUAL>::gapVec()
{
  return _phys_points_master[_qp] - _phys_points_slave[_qp];
}

template <>
DualRealVectorValue
TangentialContactLMTest<JACOBIAN>::gapVec()
{
  DualRealVectorValue gap_vec = _phys_points_master[_qp] - _phys_points_slave[_qp];
  // Here we're assuming that the user provided the x-component as the slave/master
  // variable!
  gap_vec(0).derivatives() = _u_master[_qp].derivatives() - _u_slave[_qp].derivatives();
  gap_vec(1).derivatives() =
      (*_master_disp_y_sln)[_qp].derivatives() - (*_slave_disp_y_sln)[_qp].derivatives();
  return gap_vec;
}
