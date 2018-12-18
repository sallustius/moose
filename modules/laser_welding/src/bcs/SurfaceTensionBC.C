//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceTensionBC.h"
#include "Assembly.h"

registerADMooseObject("LaserWeldingApp", SurfaceTensionBC);

defineADValidParams(SurfaceTensionBC,
                    ADIntegratedBC,
                    params.addClassDescription("Surface tension stresses.");
                    params.addRequiredParam<unsigned>("component", "The velocity component");
                    params.addParam<MaterialPropertyName>("surface_tension_name",
                                                          "surface_tension",
                                                          "The surface tension");
                    params.addParam<MaterialPropertyName>("grad_surface_tension_name",
                                                          "grad_surface_tension",
                                                          "The gradient of the surface tension"););

template <ComputeStage compute_stage>
SurfaceTensionBC<compute_stage>::SurfaceTensionBC(const InputParameters & parameters)
  : ADIntegratedBC<compute_stage>(parameters),
    _curvatures(this->_assembly.template adCurvatures<compute_stage>()),
    _component(adGetParam<unsigned>("component")),
    _surface_tension(adGetADMaterialProperty<Real>("surface_tension_name")),
    _grad_surface_tension(adGetADMaterialProperty<RealVectorValue>("grad_surface_tension_name"))
{
}

template <ComputeStage compute_stage>
ADResidual
SurfaceTensionBC<compute_stage>::computeQpResidual()
{
  return -_test[_i][_qp] *
         (2. * _curvatures[_qp] * _surface_tension[_qp] * _normals[_qp](_component) +
          _grad_surface_tension[_qp](_component) -
          _normals[_qp](_component) * (_normals[_qp] * _grad_surface_tension[_qp]));
}
