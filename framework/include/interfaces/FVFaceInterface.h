//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseMesh.h"
#include "SubProblem.h"

#include "libmesh/elem.h"

class FVFaceInterface
{
public:
  /// This codifies a set of available ways to interpolate with elem+neighbor
  /// solution information to calculate values (e.g. solution, material
  /// properties, etc.) at the face (centroid).  These methods are used in the
  /// class's interpolate functions.  Some interpolation methods are only meant
  /// to be used with advective terms (e.g. upwind), others are more
  /// generically applicable.
  enum class InterpMethod
  {
    /// (elem+neighbor)/2
    Average,
    /// weighted
    Upwind,
    // Rhie-Chow
    RhieChow
  };

  FVFaceInterface(const MooseObject * moose_object);

protected:
  /// Computes the average inteprolation
  template <typename T, typename T2, typename T3>
  void averageInterpolation(T & result, const T2 & elem, const T3 & neighbor) const;

  /// Provides interpolation of face values for non-advection-specific purposes
  /// (although it can/will still be used by advective kernels sometimes).  The
  /// interpolated value is stored in result.  This should be called when a
  /// face value needs to be computed using elem and neighbor information (e.g. a
  /// material property, solution value, etc.).  elem and neighbor represent the
  /// property/value to compute the face value for.
  template <typename T, typename T2, typename T3>
  void interpolate(InterpMethod m, T & result, const T2 & elem, const T3 & neighbor) const;

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
                   const Vector & advector) const;

  /// Return the normal at the face
  virtual const ADRealVectorValue & normal() const = 0;

  /// Return the face information
  virtual const FaceInfo & faceInfo() const = 0;

  /// Calculates and returns "grad_u dot normal" on the face to be used for
  /// diffusive terms.  If using any cross-diffusion corrections, etc. all
  /// those calculations will be handled for appropriately by this function.
  template <typename T, typename T2>
  ADReal gradUDotNormal(const T & elem_value, const T2 & neighbor_value, const FaceInfo &) const;

private:
  const SubProblem & _fvfi_subproblem;
};

template <typename T, typename T2>
ADReal
FVFaceInterface::gradUDotNormal(const T & elem_value,
                                const T2 & neighbor_value,
                                const FaceInfo & face_info) const
{
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
}

template <typename T, typename T2, typename T3, typename Vector>
void
FVFaceInterface::interpolate(
    InterpMethod m, T & result, const T2 & elem, const T3 & neighbor, const Vector & advector) const
{
  switch (m)
  {
    case InterpMethod::Average:
    {
      averageInterpolation(result, elem, neighbor);
      break;
    }

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

template <typename T, typename T2, typename T3>
void
FVFaceInterface::interpolate(InterpMethod m, T & result, const T2 & elem, const T3 & neighbor) const
{
  switch (m)
  {
    case InterpMethod::Average:
    {
      averageInterpolation(result, elem, neighbor);
      break;
    }
    default:
      mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
  }
}

template <typename T, typename T2, typename T3>
void
FVFaceInterface::averageInterpolation(T & result, const T2 & elem, const T3 & neighbor) const
{
  const auto & fi = faceInfo();
  Real elem_weight = 0, neighbor_weight = 0;

  mooseAssert(fi.neighborPtr() ? (_fvfi_subproblem.getCoordSystem(fi.elem().subdomain_id()) ==
                                  _fvfi_subproblem.getCoordSystem(fi.neighborPtr()->subdomain_id()))
                               : true,
              "Element and neighbor have different coordinate systems!");

  switch (_fvfi_subproblem.getCoordSystem(fi.elem().subdomain_id()))
  {
    case Moose::COORD_XYZ:
      elem_weight = 1., neighbor_weight = 1.;
      break;

    case Moose::COORD_RZ:
    {
      auto rz_radial_coord = _fvfi_subproblem.getAxisymmetricRadialCoord();
      elem_weight = fi.elemCentroid()(rz_radial_coord);
      neighbor_weight = fi.neighborCentroid()(rz_radial_coord);
      break;
    }

    case Moose::COORD_RSPHERICAL:
    {
      elem_weight = fi.elemCentroid()(0) * fi.elemCentroid()(0);
      neighbor_weight = fi.neighborCentroid()(0) * fi.neighborCentroid()(0);
      break;
    }

    default:
      mooseError("Unknown coordinate system");
  }

  mooseAssert(elem_weight > 0 && neighbor_weight > 0,
              "Invalid weights in FVFaceInterface::interpolate");

  result = (elem * elem_weight + neighbor * neighbor_weight) / (elem_weight + neighbor_weight);
}
