//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressDivergenceRSphericalTensors.h"
#include "ElasticityTensorTools.h"
#include "FEProblem.h"
#include "MooseMesh.h"

registerMooseObject("TensorMechanicsApp", StressDivergenceRSphericalTensors);

template <>
InputParameters
validParams<StressDivergenceRSphericalTensors>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addClassDescription(
      "Calculate stress divergence for an spherically symmetric 1D problem in polar coordinates.");
  params.addRequiredParam<unsigned int>(
      "component",
      "An integer corresponding to the direction the variable this kernel acts in. (0 "
      "for x, 1 for y, 2 for z; note in this kernel disp_x refers to the radial "
      "displacement and disp_y refers to the axial displacement.)");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceRSphericalTensors::StressDivergenceRSphericalTensors(
    const InputParameters & parameters)
  : StressDivergenceTensors(parameters)
{
  if (_component != 0)
    mooseError("Invalid component for this 1D RSpherical problem.");
}

void
StressDivergenceRSphericalTensors::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_RSPHERICAL)
    mooseError("The coordinate system in the Problem block must be set to RSPHERICAL for 1D "
               "spherically symmetric geometries.");
}

ADReal
StressDivergenceRSphericalTensors::computeQpResidual()
{
  return _grad_test[_i][_qp](0) * _stress[_qp](0, 0) +               // stress_{rr} part 1
         +(_test[_i][_qp] / _q_point[_qp](0)) * _stress[_qp](1, 1) + // stress_{\theta \theta}
         +(_test[_i][_qp] / _q_point[_qp](0)) * _stress[_qp](2, 2);  // stress_{\phi \phi}
}
