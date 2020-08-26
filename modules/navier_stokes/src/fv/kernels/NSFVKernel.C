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
#include "FVDirichletBC.h"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/vector_value.h"

registerMooseObject("NavierStokesApp", NSFVKernel);

namespace
{
ADReal
coeffCalculator(const Elem * const elem, void * context)
{
  auto * nsfv_kernel = static_cast<NSFVKernel *>(context);

  return nsfv_kernel->coeffCalculator(elem);
}
}

InputParameters
NSFVKernel::validParams()
{
  InputParameters params = FVMatAdvection::validParams();
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
NSFVKernel::coeffCalculator(const Elem * const elem)
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
    FVFluxKernel::interpolate(InterpMethod::Average, interp_v, elem_velocity, neighbor_velocity);

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
NSFVKernel::interpolate(InterpMethod m,
                        ADRealVectorValue & v,
                        const ADRealVectorValue & elem_v,
                        const ADRealVectorValue & neighbor_v)
{
  FVFluxKernel::interpolate(InterpMethod::Average, v, elem_v, neighbor_v);

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
      FVFluxKernel::interpolate(InterpMethod::Average, face_a, elem_a, neighbor_a);
      FVFluxKernel::interpolate(InterpMethod::Average,
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

ADReal
NSFVKernel::computeQpResidual()
{
  ADRealVectorValue v;
  ADReal u_interface;

  interpolate(_velocity_interp_method, v, _vel_elem[_qp], _vel_neighbor[_qp]);
  FVFluxKernel::interpolate(
      _advected_interp_method, u_interface, _adv_quant_elem[_qp], _adv_quant_neighbor[_qp], v);
  return _normal * v * u_interface;
}

#endif
