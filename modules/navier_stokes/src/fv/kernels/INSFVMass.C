//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMass.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVMass);

InputParameters
INSFVMass::validParams()
{
  InputParameters params = FVMatAdvection::validParams();
  params.addClassDescription("Enforces the divergence free condition of the velocity field.");
  params.set<MaterialPropertyName>("advected_quantity") = std::to_string(1);
  params.set<MaterialPropertyName>("vel") = NS::velocity;
  return params;
}

INSFVMass::INSFVMass(const InputParameters & params) : FVMatAdvection(params) {}
