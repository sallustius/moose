//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressDivergenceRZTensors.h"
#include "Assembly.h"
#include "ElasticityTensorTools.h"
#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", StressDivergenceRZTensors);

template <>
InputParameters
validParams<StressDivergenceRZTensors>()
{
  InputParameters params = validParams<StressDivergenceTensors>();
  params.addClassDescription(
      "Calculate stress divergence for an axisymmetric problem in cylinderical coordinates.");
  params.addRequiredParam<unsigned int>(
      "component",
      "An integer corresponding to the direction the variable this kernel acts in. (0 "
      "for x, 1 for y, 2 for z; note in this kernel disp_x refers to the radial "
      "displacement and disp_y refers to the axial displacement.)");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceRZTensors::StressDivergenceRZTensors(const InputParameters & parameters)
  : StressDivergenceTensors(parameters)
{
}

void
StressDivergenceRZTensors::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_RZ)
    mooseError("The coordinate system in the Problem block must be set to RZ for axisymmetric "
               "geometries.");
}

ADReal
StressDivergenceRZTensors::computeQpResidual()
{
  Real div = 0.0;
  if (_component == 0)
  {
    div = _grad_test[_i][_qp](0) * _stress[_qp](0, 0) +
          +(_test[_i][_qp] / _q_point[_qp](0)) * _stress[_qp](2, 2) +
          +_grad_test[_i][_qp](1) * _stress[_qp](0, 1); // stress_{rz}

    // volumetric locking correction
    if (_volumetric_locking_correction)
      div += (_avg_grad_test[_i][0] - _grad_test[_i][_qp](0) - _test[_i][_qp] / _q_point[_qp](0)) *
             (_stress[_qp].tr()) / 3.0;
  }
  else if (_component == 1)
  {
    div = _grad_test[_i][_qp](1) * _stress[_qp](1, 1) +
          +_grad_test[_i][_qp](0) * _stress[_qp](1, 0); // stress_{zr}

    // volumetric locking correction
    if (_volumetric_locking_correction)
      div += (_avg_grad_test[_i][1] - _grad_test[_i][_qp](1)) * (_stress[_qp].tr()) / 3.0;
  }
  else
    mooseError("Invalid component for this AxisymmetricRZ problem.");

  return div;
}

void
StressDivergenceRZTensors::computeAverageGradientTest()
{
  // calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(2);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      if (_component == 0)
        _avg_grad_test[_i][_component] +=
            (_grad_test[_i][_qp](_component) + _test[_i][_qp] / _q_point[_qp](0)) * _JxW[_qp] *
            _coord[_qp];
      else
        _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];
    }
    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}
