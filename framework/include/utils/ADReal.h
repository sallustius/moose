//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef ADREAL_H
#define ADREAL_H

#include "libmesh/libmesh_common.h"
#include "libmesh/compare_types.h"

namespace MetaPhysicL
{
template <typename, typename>
class DualNumber;
template <std::size_t, typename>
class NumberArray;
}

using libMesh::Real;
using MetaPhysicL::DualNumber;
using MetaPhysicL::NumberArray;

#define AD_MAX_DOFS_PER_ELEM 243

typedef DualNumber<Real, NumberArray<AD_MAX_DOFS_PER_ELEM, Real>> ADReal;

#endif
