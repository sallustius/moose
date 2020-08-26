//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NSFVFunctionBC.h"

#ifdef MOOSE_GLOBAL_AD_INDEXING

#include "MooseVariableFieldBase.h"
#include "SystemBase.h"
#include "ADReal.h"    // Moose::derivInsert
#include "MooseMesh.h" // FaceInfo methods
#include "FVDirichletBC.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/vector_value.h"

registerMooseObject("NavierStokesApp", NSFVFunctionBC);

namespace
{
ADReal
coeffCalculator(const Elem * const elem, void * context)
{
  auto * nsfv_bc = static_cast<NSFVFunctionBC *>(context);

  return nsfv_bc->coeffCalculator(elem);
}
}

InputParameters
NSFVFunctionBC::validParams()
{
  InputParameters params = FVMatAdvectionFunctionBC::validParams();
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

  params.addRequiredParam<FunctionName>("pressure_exact_solution",
                                        "The function describing the pressure exact solution.");
  return params;
}

NSFVFunctionBC::NSFVFunctionBC(const InputParameters & params)
  : FVMatAdvectionFunctionBC(params),
    _p_var(dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("pressure", 0))),
    _u_var(dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("u", 0))),
    _v_var(isParamValid("v") ? dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("v", 0))
                             : nullptr),
    _w_var(isParamValid("w") ? dynamic_cast<const MooseVariableFV<Real> *>(getFieldVar("w", 0))
                             : nullptr),
    _mu(getParam<Real>("mu")),
    _rho(getParam<Real>("rho")),
    _pressure_exact_solution(getFunction("pressure_exact_solution"))
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

ADReal
NSFVFunctionBC::coeffCalculator(const Elem * const elem)
{
  ADReal coeff = 0;

  ADRealVectorValue elem_velocity(_u_var->getElemValue(elem));

  if (_v_var)
    elem_velocity(1) = _v_var->getElemValue(elem);
  if (_w_var)
    elem_velocity(2) = _w_var->getElemValue(elem);

  for (const auto side : elem->side_index_range())
  {
    const Elem * const neighbor = elem->neighbor_ptr(side);

    bool elem_has_info = neighbor ? elem->id() < neighbor->id() : true;

    const FaceInfo * const fi = elem_has_info
                                    ? _mesh.faceInfo(elem, side)
                                    : _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));

    mooseAssert(fi, "We should have found a FaceInfo");

    const Point elem_normal = elem_has_info ? fi->normal() : Point(-fi->normal());
    const Point surface_vector = elem_normal * fi->faceArea();

    auto neighbor_value_functor = [&](const MooseVariableFV<Real> & var, unsigned int component) {
      if (neighbor)
        return var.getElemValue(neighbor);
      else
      {
        // If we don't have a neighbor, then we're along a boundary, and we may have a DirichletBC
        std::vector<FVDirichletBC *> bcs;

        _subproblem.getMooseApp()
            .theWarehouse()
            .query()
            .template condition<AttribSystem>("FVDirichletBC")
            .template condition<AttribThread>(_tid)
            .template condition<AttribBoundaries>(_face_info->boundaryIDs())
            .template condition<AttribVar>(var.number())
            .template condition<AttribSysNum>(var.sys().number())
            .queryInto(bcs);
        mooseAssert(bcs.size() <= 1, "cannot have multiple dirichlet BCs on the same boundary");

        bool has_dirichlet_bc = bcs.size() > 0;

        if (has_dirichlet_bc)
        {
          const FVDirichletBC & bc = *bcs[0];

          // Linear interpolation: face_value = (elem_value + neighbor_value) / 2
          return 2. * bc.boundaryValue(*_face_info) - elem_velocity(component);
        }
        else
          // No DirichletBC so we'll implicitly apply a zero gradient condition and assume that the
          // face value is equivalent to the element value
          return elem_velocity(component);
      }
    };

    ADRealVectorValue neighbor_velocity(neighbor_value_functor(*_u_var, 0));
    if (_v_var)
      elem_velocity(1) = neighbor_value_functor(*_v_var, 1);
    if (_w_var)
      elem_velocity(2) = neighbor_value_functor(*_w_var, 2);

    ADRealVectorValue interp_v;
    FVFluxBC::interpolate(InterpMethod::Average, interp_v, elem_velocity, neighbor_velocity);

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

void
NSFVFunctionBC::interpolate(InterpMethod m,
                            ADRealVectorValue & v_face,
                            const ADRealVectorValue & elem_v,
                            const RealVectorValue & ghost_v)
{
  FVFluxBC::interpolate(InterpMethod::Average, v_face, elem_v, ghost_v);

  if (m == InterpMethod::RhieChow)
  {
    // Get pressure gradient for the elem
    const VectorValue<ADReal> & grad_p_elem = _p_var->adGradSln(&_face_info->elem());

    // Get pressure gradient for the ghost
    RealVectorValue grad_p_ghost = _pressure_exact_solution.gradient(
        _t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid());

    // Uncorrected face pressure gradient
    auto unc_grad_p = (grad_p_elem + grad_p_ghost) / 2.;

    // Now perform correction

    // Get pressure value for the elem
    ADReal p_elem = _p_var->getElemValue(&_face_info->elem());

    // Get pressure value for the ghost
    Real p_ghost = _pressure_exact_solution.value(
        _t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid());

    auto d_cf = 2. * (_face_info->faceCentroid() - _face_info->elemCentroid());
    auto d_cf_norm = d_cf.norm();
    auto e_cf = d_cf / d_cf_norm;

    // Corrected face pressure gradient
    auto grad_p = unc_grad_p + ((p_ghost - p_elem) / d_cf_norm - unc_grad_p * e_cf) * e_cf;

    // Now we need to perform the computations of D
    // I don't know how I would want to do computation of the a coefficient on a ghost cell. I would
    // have to essentially create an entire fictional element with defined geometric locations of
    // the faces in order to compute inward advective flux and diffusive flux. For now I'm going to
    // try not doing that and just use the a coeff of the elem
    const ADReal & face_a = _p_var->adCoeff(&_face_info->elem(), this, &::coeffCalculator);
    Real face_volume = _face_info->elemVolume();

    mooseAssert(face_a > 0, "face_a should be greater than zero");
    const ADReal face_D = face_volume / face_a;

    // perform the pressure correction
    v_face -= face_D * (grad_p - unc_grad_p);
  }
}

ADReal
NSFVFunctionBC::computeQpResidual()
{
  ADReal flux_var_face;
  ADRealVectorValue v_face;

  Real flux_var_ghost = _flux_variable_exact_solution.value(
      _t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid());
  RealVectorValue v_ghost(
      _vel_x_exact_solution.value(_t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid()),
      _vel_y_exact_solution ? _vel_y_exact_solution->value(
                                  _t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid())
                            : 0,
      _vel_z_exact_solution ? _vel_z_exact_solution->value(
                                  _t, 2. * _face_info->faceCentroid() - _face_info->elemCentroid())
                            : 0);

  interpolate(_velocity_interp_method, v_face, _vel[_qp], v_ghost);
  FVFluxBC::interpolate(
      _advected_interp_method, flux_var_face, _adv_quant[_qp], flux_var_ghost, v_face);
  return _normal * v_face * flux_var_face;
}

#endif
