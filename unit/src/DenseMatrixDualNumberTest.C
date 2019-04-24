//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "DualReal.h"
#include "MaterialProperty.h"
#include "libmesh/dense_matrix.h"
#include "metaphysicl/dualnumberarray.h"
#include <vector>

using libMesh::DenseMatrix;

TEST(DenseMatrixDualNumber, DenseMatrixDualNumber)
{
  MaterialProperty<std::vector<DenseMatrix<Real>>> * property =
      new ADMaterialPropertyObject<std::vector<DenseMatrix<Real>>>(true);

  property->resize(1);
  (*property)[0].resize(1);

  delete property;
}
