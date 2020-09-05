//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVFunctionFluxBC.h"

registerMooseObject("MooseApp", FVFunctionFluxBC);

InputParameters
FVFunctionFluxBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addRequiredParam<FunctionName>("flux_function",
                                        "The value of the flux crossing the boundary.");
  return params;
}

FVFunctionFluxBC::FVFunctionFluxBC(const InputParameters & parameters)
  : FVFluxBC(parameters), _function(getFunction("flux_function"))
{
}

ADReal
FVFunctionFluxBC::computeQpResidual()
{
  return _flux_function.value(_t, _face_info->faceCentroid());
}
