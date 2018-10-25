//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GaussianWeldEnergyFluxBC.h"

registerMooseObject("LaserWeldingApp", GaussianWeldEnergyFluxBC);

template <>
InputParameters
validParams<GaussianWeldEnergyFluxBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("reff",
                                "The effective radius describing the radial distribution of the "
                                "beam energy. This should be non-dimensional.");
  params.addRequiredParam<Real>("F0", "The average heat flux of the laser");
  params.addRequiredParam<Real>("R", "The beam radius");
  params.addRequiredParam<RealVectorValue>(
      "beam_coords", "The coordinates of the center of the beam at the surface contact point");
  return params;
}

GaussianWeldEnergyFluxBC::GaussianWeldEnergyFluxBC(const InputParameters & params)
  : IntegratedBC(params),
    _reff(getParam<Real>("reff")),
    _F0(getParam<Real>("F0")),
    _R(getParam<Real>("R")),
    _beam_coords(getParam<RealVectorValue>("beam_coords"))
{
}

Real
GaussianWeldEnergyFluxBC::computeQpResidual()
{
  auto r = (_q_point[_qp] - _beam_coords).norm();
  return -_test[_i][_qp] * 2. * _reff * _F0 * std::exp(-_reff * r * r / (_R * _R));
}
