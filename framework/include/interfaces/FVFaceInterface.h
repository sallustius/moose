//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVUtils.h"
#include "MooseVariableFV.h"

class FVFaceInterface
{
public:
  using InterpMethod = Moose::InterpMethod;

  FVFaceInterface() = default;

  /// Provides interpolation of face values for non-advection-specific purposes
  /// (although it can/will still be used by advective kernels sometimes).  The
  /// interpolated value is stored in result.  This should be called when a
  /// face value needs to be computed using elem and neighbor information (e.g. a
  /// material property, solution value, etc.).  elem and neighbor represent the
  /// property/value to compute the face value for.
  template <typename T, typename T2, typename T3>
  void interpolate(InterpMethod m, T & result, const T2 & elem, const T3 & neighbor) const
  {
    switch (m)
    {
      case InterpMethod::Average:
        result = (elem + neighbor) * 0.5;
        break;
      default:
        mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
    }
  }

  /// Provides interpolation of face values for advective flux kernels.  This
  /// should be called by advective kernels when a u_face value is needed from
  /// u_elem and u_neighbor.  The interpolated value is stored in result.  elem
  /// and neighbor represent the property/value being advected in the elem and
  /// neighbor elements respectively.  advector represents the vector quantity at
  /// the face that is doing the advecting (e.g. the flow velocity at the
  /// face); this value often will have been computed using a call to the
  /// non-advective interpolate function.
  template <typename T, typename T2, typename T3, typename Vector>
  void interpolate(InterpMethod m,
                   T & result,
                   const T2 & elem,
                   const T3 & neighbor,
                   const Vector & advector) const
  {
    switch (m)
    {
      case InterpMethod::Average:
        result = (elem + neighbor) * 0.5;
        break;
      case InterpMethod::Upwind:
        if (advector * normal() > 0)
          result = elem;
        else
          result = neighbor;
        break;
      default:
        mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
    }
  }

  /// Return the normal at the face
  virtual const ADRealVectorValue & normal() const = 0;

  /// Calculates and returns "grad_u dot normal" on the face to be used for
  /// diffusive terms.  If using any cross-diffusion corrections, etc. all
  /// those calculations will be handled for appropriately by this function.
  template <typename T, typename T2>
  ADReal gradUDotNormal(const T & elem_value,
                        const T2 & neighbor_value,
                        const FaceInfo & fi,
                        const MooseVariableFV<Real> & fv_var) const;
};

template <typename T, typename T2>
ADReal
FVFaceInterface::gradUDotNormal(
#ifdef MOOSE_GLOBAL_AD_INDEXING
    const T &, const T2 &, const FaceInfo & face_info, const MooseVariableFV<Real> & fv_var
#else
    const T & elem_value,
    const T2 & neighbor_value,
    const FaceInfo & face_info,
    const MooseVariableFV<Real> &
#endif
) const
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  return fv_var.adGradSln(face_info) * normal();

#else
  // We compute "grad_u dot _normal" by assuming the mesh is orthogonal, and
  // recognizing that it is equivalent to delta u between the two cell
  // centroids but for one unit in the normal direction.  We know delta u for
  // the length between cell centroids (u_neighbor - u_elem) and then we just
  // divide that by the distance between the centroids to convert it to delta
  // u for one unit in the normal direction.  Because the _normal vector is
  // defined to be outward from the elem element, u_neighbor-u_elem gives delta u
  // when moving in the positive normal direction.  So we divide by the
  // (positive) distance between centroids because one unit in the normal
  // direction is always positive movement.
  ADReal dudn =
      (neighbor_value - elem_value) /
      (face_info.neighborPtr() ? face_info.neighborCentroid() - face_info.elemCentroid()
                               : 2. * (face_info.faceCentroid() - face_info.elemCentroid()))
          .norm();

  // TODO: need to apply cross-diffusion correction factor here.  This
  // currently is only correct if the vector between the elem-neighbor cell
  // centroids is parallel to the normal vector.
  return dudn;
#endif
}
