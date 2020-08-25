//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVMatAdvection.h"
#include "FVMatAdvectionFunctionBC.h"

#ifdef MOOSE_GLOBAL_AD_INDEXING

template <typename>
class MooseVariableFV;

namespace libMesh
{
class DofMap;
template <typename>
class NumericVector;
template <typename>
class VectorValue;
}

template <typename T>
class NSFVBase : public T
{
public:
  static InputParameters validParams();
  NSFVBase(const InputParameters & params);

  ADReal coeffCalculator(const Elem * const elem);

protected:
  /**
   * interpolation overload for the velocity
   */
  void interpolate(InterpMethod m,
                   ADRealVectorValue & interp_v,
                   const ADRealVectorValue & elem_v,
                   const ADRealVectorValue & neighbor_v);

  ADReal computeQpResidual() override;

private:
  /// pressure variable
  const MooseVariableFV<Real> * const _p_var;
  /// x-velocity
  const MooseVariableFV<Real> * const _u_var;
  /// y-velocity
  const MooseVariableFV<Real> * const _v_var;
  /// z-velocity
  const MooseVariableFV<Real> * const _w_var;

  /// The viscosity
  const Real _mu;
  /// The density
  const Real _rho;

  /// The interpolation method to use for the velocity
  InterpMethod _velocity_interp_method;
};

template <typename T>
InputParameters
NSFVBase<T>::validParams()
{
  InputParameters params = T::validParams();
  params.addRequiredCoupledVar("pressure", "The pressure variable.");
  params.addRequiredCoupledVar("u", "The velocity in the x direction.");
  params.addCoupledVar("v", "The velocity in the y direction.");
  params.addCoupledVar("w", "The velocity in the z direction.");

  MooseEnum velocity_interp_method("average rc", "rc");

  params.addParam<MooseEnum>(
      "velocity_interp_method",
      velocity_interp_method,
      "The interpolation to use for the velocity. Options are "
      "'average' and 'rc' which stands for Rhie-Chow. The default is Rhie-Chow.");

  params.addParam<Real>("mu", 1, "The viscosity");
  params.addParam<Real>("rho", 1, "The density");
  return params;
}

template <typename T>
NSFVBase<T>::NSFVBase(const InputParameters & params)
  : T(params),
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

  const auto & velocity_interp_method = getParam<MooseEnum>("velocity_interp_method");
  if (velocity_interp_method == "average")
    _velocity_interp_method = InterpMethod::Average;
  else if (velocity_interp_method == "rc")
    _velocity_interp_method = InterpMethod::RhieChow;
  else
    mooseError("Unrecognized interpolation type ",
               static_cast<std::string>(velocity_interp_method));
}

template <typename T>
ADReal
NSFVBase<T>::coeffCalculator(const Elem * const elem)
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
    T::interpolate(InterpMethod::Average, interp_v, elem_velocity, neighbor_velocity);

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

template <typename T>
void
NSFVBase<T>::interpolate(InterpMethod m,
                         ADRealVectorValue & v,
                         const ADRealVectorValue & elem_v,
                         const ADRealVectorValue & neighbor_v)
{
  T::interpolate(InterpMethod::Average, v, elem_v, neighbor_v);

  if (m == InterpMethod::RhieChow)
  {
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
      T::interpolate(InterpMethod::Average, face_a, elem_a, neighbor_a);
      T::interpolate(InterpMethod::Average,
                     face_volume,
                     _face_info->elemVolume(),
                     _face_info->neighborVolume());
    }

    mooseAssert(face_a > 0, "face_a should be greater than zero");
    const ADReal face_D = face_volume / face_a;

    // perform the pressure correction
    v -= face_D * (grad_p - unc_grad_p);
  }
}

template <typename T>
ADReal
NSFVBase<T>::computeQpResidual()
{
  ADRealVectorValue v;
  ADReal u_interface;

  interpolate(_velocity_interp_method, v, _vel_elem[_qp], _vel_neighbor[_qp]);
  T::interpolate(
      _advected_interp_method, u_interface, _adv_quant_elem[_qp], _adv_quant_neighbor[_qp], v);
  return _normal * v * u_interface;
}

#endif
