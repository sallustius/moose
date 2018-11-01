//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledPenaltyInterfaceDiffusion.h"

registerMooseObject("MooseApp", CoupledPenaltyInterfaceDiffusion);

template <>
InputParameters
validParams<CoupledPenaltyInterfaceDiffusion>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addRequiredParam<Real>(
      "penalty", "The penalty that penalizes jump between master and neighbor variables.");
  params.addCoupledVar("master_coupled_var", "The coupled variable on the master side");
  params.addCoupledVar("slave_coupled_var", "The coupled variable on the slave side");
  return params;
}

CoupledPenaltyInterfaceDiffusion::CoupledPenaltyInterfaceDiffusion(
    const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _penalty(getParam<Real>("penalty")),
    _master_coupled_value(isCoupled("master_coupled_var") ? coupledValue("master_coupled_var")
                                                          : _var.sln()),
    _slave_coupled_value(isCoupled("slave_coupled_var") ? coupledNeighborValue("slave_coupled_var")
                                                        : _neighbor_var.slnNeighbor())
{
}

Real
CoupledPenaltyInterfaceDiffusion::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * _penalty * (_master_coupled_value[_qp] - _slave_coupled_value[_qp]);
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] * -_penalty *
          (_master_coupled_value[_qp] - _slave_coupled_value[_qp]);
      break;
  }

  return r;
}
