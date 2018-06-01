//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LM.h"
#include "PenetrationInfo.h"
#include "PenetrationLocator.h"

#include "libmesh/string_to_enum.h"

registerMooseObject("ContactApp", LM);

template <>
InputParameters
validParams<LM>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addRequiredParam<BoundaryName>("slave", "The boundary ID associated with the slave side");
  params.addRequiredParam<BoundaryName>("master",
                                        "The boundary ID associated with the master side");
  MooseEnum orders("FIRST SECOND THIRD FOURTH", "FIRST");
  params.addParam<MooseEnum>("order", orders, "The finite element order used for projections");
  params.set<bool>("use_displaced_mesh") = true;
  params.addCoupledVar("disp_x", "The x displacement");
  params.addCoupledVar("disp_y", "The y displacement");
  params.addCoupledVar("disp_z", "The z displacement");
  return params;
}

LM::LM(const InputParameters & parameters)
  : NodalKernel(parameters),
    _penetration_locator(
        getPenetrationLocator(getParam<BoundaryName>("master"),
                              getParam<BoundaryName>("slave"),
                              Utility::string_to_enum<Order>(getParam<MooseEnum>("order")))),
    _disp_x_id(coupled("disp_x")),
    _disp_y_id(coupled("disp_y")),
    _disp_z_id(coupled("disp_z"))
{
}

Real
LM::computeQpResidual()
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo != NULL)
    {
      Real a = -pinfo->_distance;
      Real b = _u[_qp];
      return a + b - std::sqrt(a * a + b * b);
    }
  }
  return 0;
}

Real
LM::computeQpJacobian()
{
  std::map<dof_id_type, PenetrationInfo *>::iterator found =
      _penetration_locator._penetration_info.find(_current_node->id());
  if (found != _penetration_locator._penetration_info.end())
  {
    PenetrationInfo * pinfo = found->second;
    if (pinfo != NULL)
    {
      Real a = pinfo->_distance;
      Real b = _u[_qp];
      return 1. + b / std::sqrt(a * a + b * b);
    }
  }
  return 0;
}

Real
LM::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _disp_x_id || jvar == _disp_y_id || jvar == _disp_z_id)
  {
    std::map<dof_id_type, PenetrationInfo *>::iterator found =
        _penetration_locator._penetration_info.find(_current_node->id());
    if (found != _penetration_locator._penetration_info.end())
    {
      PenetrationInfo * pinfo = found->second;
      if (pinfo != NULL)
      {
        const Elem * master_elem = pinfo->_elem;
        unsigned master_side = pinfo->_side_num;

        /*
         * Make sure proper values of test and shape functions are in _test and _phi
         */
        std::vector<Point> points = {pinfo->_closest_point};

        // reinit variables on the master element's face at the contact point
        _fe_problem.setNeighborSubdomainID(master_elem, 0);
        _fe_problem.reinitNeighborPhys(master_elem, master_side, points, 0);

        MooseVariableFEBase & var = _sys.getVariable(0, jvar);
        std::vector<dof_id_type> master_dof_indices;
        var.getDofIndices(master_elem, dof_indices);
      }
    }
  }
}
