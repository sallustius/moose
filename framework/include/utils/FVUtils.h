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
#include "MooseError.h"
#include "SubProblem.h"
#include "Assembly.h"
#include "libmesh/elem.h"

namespace Moose
{
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

template <typename ActionFunctor>
void
loopOverElemFaceInfo(const Elem & elem,
                     const MooseMesh & mesh,
                     const SubProblem & subproblem,
                     ActionFunctor & act)
{
  for (const auto side : elem.side_index_range())
  {
    const Elem * const neighbor = elem.neighbor_ptr(side);

    bool elem_has_info = neighbor ? (elem.id() < neighbor->id()) : true;

    const FaceInfo * const fi = elem_has_info
                                    ? mesh.faceInfo(&elem, side)
                                    : mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(&elem));

    mooseAssert(fi, "We should have found a FaceInfo");

    const Point elem_normal = elem_has_info ? fi->normal() : Point(-fi->normal());

    mooseAssert(neighbor ? subproblem.getCoordSystem(elem.subdomain_id()) ==
                               subproblem.getCoordSystem(neighbor->subdomain_id())
                         : true,
                "Coordinate systems must be the same between element and neighbor");

    Real coord;
    coordTransformFactor(subproblem, elem.subdomain_id(), fi->faceCentroid(), coord);

    const Point surface_vector = elem_normal * fi->faceArea() * coord;

    act(elem, neighbor, fi, surface_vector, coord, elem_has_info);
  }
}
}
