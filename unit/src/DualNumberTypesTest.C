//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "MooseTypes.h"
#include <iostream>

const double tol = 1e-8;

TEST(DualNumberTypesTest, sample)
{
  DualNumber<Real, NumberArray<2, Real>> scalar_ad_prop;
  DualNumber<VectorValue<Real>, NumberArray<2, VectorValue<Real>>> vector_ad_prop;
  DualNumber<TensorValue<Real>, NumberArray<2, TensorValue<Real>>> tensor_ad_prop;

  Real scalar_reg_prop;
  RealVectorValue vector_reg_prop;
  RealTensorValue tensor_reg_prop;

  scalar_ad_prop = 2.;
  scalar_ad_prop.derivatives()[0] = 1.;
  scalar_ad_prop.derivatives()[1] = 2.;
  scalar_ad_prop *= scalar_ad_prop;
  EXPECT_NEAR(scalar_ad_prop.value(), 4., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[0], 4., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[1], 8., tol);

  scalar_ad_prop = scalar_ad_prop * scalar_ad_prop;
  EXPECT_NEAR(scalar_ad_prop.value(), 16., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[0], 32., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[1], 64., tol);

  vector_ad_prop = scalar_ad_prop * RealVectorValue(1., 1., 1.);
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    EXPECT_NEAR(vector_ad_prop.value()(i), 16., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[0](i), 32., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[1](i), 64., tol);
  }
  scalar_ad_prop = vector_ad_prop * vector_ad_prop;
  EXPECT_NEAR(scalar_ad_prop.value(), 768., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[0], 3072., tol);
  EXPECT_NEAR(scalar_ad_prop.derivatives()[1], 6144., tol);

  vector_ad_prop *= 3.;
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    EXPECT_NEAR(vector_ad_prop.value()(i), 48., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[0](i), 96., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[1](i), 192., tol);
  }
  vector_ad_prop = 4. * vector_ad_prop + vector_ad_prop * 5.;
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    EXPECT_NEAR(vector_ad_prop.value()(i), 432., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[0](i), 864., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[1](i), 1728., tol);
  }

  tensor_ad_prop = scalar_ad_prop * RealTensorValue(1., 1., 1., 1., 1., 1., 1., 1., 1.);
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    for (decltype(LIBMESH_DIM) j = 0; j != LIBMESH_DIM; ++j)
    {
      EXPECT_NEAR(tensor_ad_prop.value()(i, j), 768., tol);
      EXPECT_NEAR(tensor_ad_prop.derivatives()[0](i, j), 3072., tol);
      EXPECT_NEAR(tensor_ad_prop.derivatives()[1](i, j), 6144., tol);
    }
  }
  tensor_ad_prop *= 1. / 768.;
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    for (decltype(LIBMESH_DIM) j = 0; j != LIBMESH_DIM; ++j)
    {
      EXPECT_NEAR(tensor_ad_prop.value()(i, j), 1., tol);
      EXPECT_NEAR(tensor_ad_prop.derivatives()[0](i, j), 4., tol);
      EXPECT_NEAR(tensor_ad_prop.derivatives()[1](i, j), 8., tol);
      tensor_ad_prop.value()(i, j) = i + j;
      tensor_ad_prop.derivatives()[0](i, j) = 2. * i + 4. * j;
      tensor_ad_prop.derivatives()[1](i, j) = 3. * i + 5. * j;
    }
  }
  vector_ad_prop *= 1. / 432.;
  for (decltype(LIBMESH_DIM) i = 0; i != LIBMESH_DIM; ++i)
  {
    EXPECT_NEAR(vector_ad_prop.value()(i), 1., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[0](i), 2., tol);
    EXPECT_NEAR(vector_ad_prop.derivatives()[1](i), 4., tol);
    vector_ad_prop.value()(i) = 10. - i;
    vector_ad_prop.derivatives()[0](i) = 10. - 2. * i;
    vector_ad_prop.derivatives()[1](i) = 10. - 3. * i;
  }

  vector_ad_prop = tensor_ad_prop * vector_ad_prop;

  tensor_ad_prop *= tensor_ad_prop;
  tensor_ad_prop = tensor_ad_prop * tensor_ad_prop;
  tensor_ad_prop = 7. * tensor_ad_prop + tensor_ad_prop * 8.;
  tensor_ad_prop = 9. * tensor_ad_prop - tensor_ad_prop * 10.;

  scalar_reg_prop = scalar_ad_prop;
  vector_reg_prop = vector_ad_prop;
  tensor_reg_prop = tensor_ad_prop;
}
