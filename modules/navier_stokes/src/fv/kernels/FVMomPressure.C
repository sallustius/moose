//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVMomPressure.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", FVMomPressure);

InputParameters
FVMomPressure::validParams()
{
  InputParameters params = FVMatAdvection::validParams();
  params.addClassDescription(
      "Introduces the coupled pressure term into the Navier-Stokes momentum equation.");
  MooseEnum momentum_component("x=0 y=1 z=2", "x");
  params.addParam<MooseEnum>("momentum_component",
                             momentum_component,
                             "The component of the momentum equation that this BC applies to.");
  params.set<MaterialPropertyName>("advected_quantity") = NS::pressure;
  params.set<MooseEnum>("advected_interp_method") = "average";
  return params;
}

FVMomPressure::FVMomPressure(const InputParameters & params)
  : FVMatAdvection(params), _index(getParam<MooseEnum>("momentum_component"))
{
}

ADReal
FVMomPressure::computeQpResidual()
{
  ADRealVectorValue v;
  ADReal p_interface;
  interpolate(InterpMethod::Average, v, _vel_elem[_qp], _vel_neighbor[_qp]);
  interpolate(
      _advected_interp_method, p_interface, _adv_quant_elem[_qp], _adv_quant_neighbor[_qp], v);

  // The difference between this residual computation and that in FVMatAdvection is:
  // _normal * v -> _normal(_index)
  return _normal(_index) * p_interface;
}
