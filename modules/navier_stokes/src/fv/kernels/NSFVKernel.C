//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NSFVKernel.h"

#ifdef MOOSE_GLOBAL_AD_INDEXING

#include "MooseVariableFieldBase.h"
#include "SystemBase.h"
#include "ADReal.h"    // Moose::derivInsert
#include "MooseMesh.h" // FaceInfo methods

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/vector_value.h"

registerMooseObject("NavierStokesApp", NSFVKernel);

ADReal
coeffCalculator(const Elem & elem, void * context)
{
  auto * nsfv_kernel = static_cast<NSFVKernel *>(context);

  return nsfv_kernel->coeffCalculator(elem);
}

InputParameters
NSFVKernel::validParams()
{
  InputParameters params = FVMatAdvection::validParams();
  params.addRequiredCoupledVar("pressure", "The pressure variable.");
  params.addRequiredCoupledVar("u", "The velocity in the x direction.");
  params.addCoupledVar("v", "The velocity in the y direction.");
  params.addCoupledVar("w", "The velocity in the z direction.");

  params.suppressParameter<MooseEnum>("advected_interp_method");

  params.addParam<Real>("mu", 1, "The viscosity");
  params.addParam<Real>("rho", 1, "The density");
  return params;
}

NSFVKernel::NSFVKernel(const InputParameters & params)
  : FVMatAdvection(params),
    _p_var(dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("pressure", 0))),
    _u_var(dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("u", 0))),
    _v_var(isParamValid("v") ? dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("v", 0))
                             : nullptr),
    _w_var(isParamValid("w") ? dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("w", 0))
                             : nullptr),
    _mu(getParam<Real>("mu")),
    _rho(getParam<Real>("rho"))
{
  if (!_p_var)
    paramError("pressure", "the pressure must be a finite volume variable.");

  if (!_u_var)
    paramError("u", "the u velocity must be a finite volume variable.");

  if (_mesh.dimension() >= 2 && !_v_var)
    mooseError("In two-dimensions, the v velocity must be supplied and it must be a finite volume "
               "variable.");

  if (_mesh.dimension() >= 3 && !isParamValid("w"))
    mooseError("In three-dimensions, the w velocity must be supplied and it must be a finite "
               "volume variable.");
}

ADReal
NSFVKernel::coeffCalculator(const Elem & elem)
{
  ADReal coeff = 0;

  ADRealVectorValue elem_velocity(_u_var->getElemValue(&elem));

  if (_v_var)
    elem_velocity(1) = _v_var->getElemValue(&elem);
  if (_w_var)
    elem_velocity(2) = _w_var->getElemValue(&elem);

  for (const auto side : elem.side_index_range())
  {
    const Elem * const neighbor = elem.neighbor_ptr(side);
    if (!neighbor)
      continue;

    bool elem_has_info = elem.id() < neighbor->id();

    const FaceInfo * const fi =
        elem_has_info ? _mesh.faceInfo(&elem, side)
                      : _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(&elem));

    mooseAssert(fi, "We should have found a FaceInfo");

    const Point elem_normal = elem_has_info ? fi->normal() : Point(-fi->normal());
    const Point surface_vector = elem_normal * fi->faceArea();

    ADRealVectorValue neighbor_velocity(_u_var->getElemValue(neighbor));
    if (_v_var)
      elem_velocity(1) = _v_var->getElemValue(neighbor);
    if (_w_var)
      elem_velocity(2) = _w_var->getElemValue(neighbor);

    ADRealVectorValue interp_v;
    interpolate(InterpMethod::Average, interp_v, elem_velocity, neighbor_velocity);

    const ADReal mass_flow = _rho * interp_v * surface_vector;

    // We are upwinding, so we only sum into the coefficient if the mass flow rate is negative,
    // indicating *inflow*
    if (mass_flow < 0)
      coeff += -mass_flow;

    // Now add the viscous flux
    coeff += _mu * fi->faceArea() / (fi->elemCentroid() - fi->neighborCentroid()).norm();
  }

  return coeff;
}

ADRealVectorValue
NSFVKernel::interpolateVelocity()
{
  ADRealVectorValue v;

  // Currently only Average is supported for the velocity
  interpolate(InterpMethod::Average, v, _vel_elem[_qp], _vel_neighbor[_qp]);

  // Get pressure gradient
  const VectorValue<ADReal> & grad_p = _p_var->adGradSln(*_face_info);

  // Get uncorrected pressure gradient
  const VectorValue<ADReal> & unc_grad_p = _p_var->uncorrectedAdGradSln(*_face_info);

  // Now we need to perform the computations of D
  const ADReal & elem_a = _p_var->adCoeff(&_face_info->elem(), this, &::coeffCalculator);
  ADReal face_a = elem_a;
  Real face_volume = _face_info->elemVolume();

  if (_face_info->neighborPtr())
  {
    const ADReal & neighbor_a =
        _p_var->adCoeff(_face_info->neighborPtr(), this, &::coeffCalculator);
    interpolate(InterpMethod::Average, face_a, elem_a, neighbor_a);
    interpolate(
        InterpMethod::Average, face_volume, _face_info->elemVolume(), _face_info->neighborVolume());
  }

  mooseAssert(face_a > 0, "face_a should be greater than zero");
  const ADReal face_D = face_volume / face_a;

  // perform the pressure correction
  v -= face_D * (grad_p - unc_grad_p);

  return v;
}

ADReal
NSFVKernel::computeQpResidual()
{
  ADRealVectorValue v = interpolateVelocity();
  ADReal u_interface;

  interpolate(InterpMethod::Upwind, u_interface, _adv_quant_elem[_qp], _adv_quant_neighbor[_qp], v);
  return _normal * v * u_interface;
}

#endif
