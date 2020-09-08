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
#include "NSFVUtils.h"

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
  params += NSFVBase::validParams();
  return params;
}

NSFVKernel::NSFVKernel(const InputParameters & params) : FVMatAdvection(params), NSFVBase(params) {}

ADReal
NSFVKernel::coeffCalculator(const Elem * const elem)
{
  return NS::coeffCalculator(elem, *this);
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

    mooseAssert(_face_info->neighborPtr()
                    ? _subproblem.getCoordSystem(_face_info->elem().subdomain_id()) ==
                          _subproblem.getCoordSystem(_face_info->neighborPtr()->subdomain_id())
                    : true,
                "Coordinate systems must be the same between element and neighbor");

    Real coord;
    coordTransformFactor(
        _subproblem, _face_info->elem().subdomain_id(), _face_info->elemCentroid(), coord);

    Real elem_volume = _face_info->elemVolume() * coord;
    Real face_volume = elem_volume;

    if (_face_info->neighborPtr())
    {
      const ADReal & neighbor_a =
          _p_var->adCoeff(_face_info->neighborPtr(), this, &::coeffCalculator);
      FVFluxKernel::interpolate(InterpMethod::Average, face_a, elem_a, neighbor_a);

      coordTransformFactor(_subproblem,
                           _face_info->neighborPtr()->subdomain_id(),
                           _face_info->neighborCentroid(),
                           coord);
      Real neighbor_volume = _face_info->neighborVolume() * coord;
      FVFluxKernel::interpolate(InterpMethod::Average, face_volume, elem_volume, neighbor_volume);
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
