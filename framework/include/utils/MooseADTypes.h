#ifndef MOOSEADTYPES_H
#define MOOSEADTYPES_H

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"

template <typename>
class MooseArray;

// The 100 here is for how many DoFs there are per element.
#define AD_MAX_DOFS_PER_ELEM 100
using MetaPhysicL::DualNumber;
using MetaPhysicL::NumberArray;

template <typename T>
using ScalarDN = DualNumber<T, NumberArray<AD_MAX_DOFS_PER_ELEM, T>>;
template <typename T>
using VectorDN = DualNumber<VectorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>>>;
template <typename T>
using TensorDN = DualNumber<TensorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>>>;

typedef ScalarDN<Real> ADReal;

template <typename T>
VectorDN<T> operator*(const VectorDN<T> & vec, const T & scalar)
{
  VectorValue<T> value = vec.value() * scalar;
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec.derivatives()[i] * scalar;
  return {value, derivatives};
}

template <typename T>
VectorDN<T> operator*(const T & scalar, const VectorDN<T> & vec)
{
  VectorValue<T> value = vec.value() * scalar;
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec.derivatives()[i] * scalar;
  return {value, derivatives};
}

template <typename T>
VectorDN<T> operator*(const VectorValue<T> & vec, const ScalarDN<T> & scalar)
{
  VectorValue<T> value = vec * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec * scalar.derivatives()[i];
  return {value, derivatives};
}

template <typename T>
VectorDN<T> operator*(const ScalarDN<T> & scalar, const VectorValue<T> & vec)
{
  VectorValue<T> value = vec * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec * scalar.derivatives()[i];
  return {value, derivatives};
}

template <typename T>
VectorDN<T> operator*(const VectorDN<T> & vec, const ScalarDN<T> & scalar)
{
  VectorValue<T> value = vec.value() * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec.value() * scalar.derivatives()[i] + vec.derivatives()[i] * scalar.value();
  return {value, derivatives};
}

template <typename T>
VectorDN<T> operator*(const ScalarDN<T> & scalar, const VectorDN<T> & vec)
{
  VectorValue<T> value = vec.value() * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = vec.value() * scalar.derivatives()[i] + vec.derivatives()[i] * scalar.value();
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const TensorDN<T> & tensor, const T & scalar)
{
  TensorValue<T> value = tensor.value() * scalar;
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = tensor.derivatives()[i] * scalar;
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const T & scalar, const TensorDN<T> & tensor)
{
  TensorValue<T> value = tensor.value() * scalar;
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = tensor.derivatives()[i] * scalar;
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const TensorValue<T> & tensor, const ScalarDN<T> & scalar)
{
  TensorValue<T> value = tensor * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = tensor * scalar.derivatives()[i];
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const ScalarDN<T> & scalar, const TensorValue<T> & tensor)
{
  TensorValue<T> value = tensor * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = tensor * scalar.derivatives()[i];
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const TensorDN<T> & tensor, const ScalarDN<T> & scalar)
{
  TensorValue<T> value = tensor.value() * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] =
        tensor.value() * scalar.derivatives()[i] + tensor.derivatives()[i] * scalar.value();
  return {value, derivatives};
}

template <typename T>
TensorDN<T> operator*(const ScalarDN<T> & scalar, const TensorDN<T> & tensor)
{
  TensorValue<T> value = tensor.value() * scalar.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] =
        tensor.value() * scalar.derivatives()[i] + tensor.derivatives()[i] * scalar.value();
  return {value, derivatives};
}

template <typename T>
ScalarDN<T> operator*(const VectorDN<T> & dn, const VectorValue<T> & vec2)
{
  const VectorValue<T> & vec1 = dn.value();
  T value = vec1 * vec2;
  NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = dn.derivatives()[i] * vec2;
  return {value, derivatives};
}

template <typename T>
ScalarDN<T> operator*(const VectorValue<T> & vec2, const VectorDN<T> & dn)
{
  const VectorValue<T> & vec1 = dn.value();
  T value = vec1 * vec2;
  NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = dn.derivatives()[i] * vec2;
  return {value, derivatives};
}

template <typename T>
ScalarDN<T> operator*(const VectorDN<T> & dn1, const VectorDN<T> & dn2)
{
  T value = dn1.value() * dn2.value();
  NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
  for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
    derivatives[i] = dn1.derivatives()[i] * dn2.value() + dn1.value() * dn2.derivatives()[i];
  return {value, derivatives};
}

// template <typename MatTemplateType, typename ScalarType, template MatType>
//     > auto operator*(const typename MatType<MatTemplateType> & mat, const ScalarType & scalar)
//           -> MatType<decltype(mat(0, 0) * scalar)>
// {
//   MatType<decltype(mat(0, 0) * scalar)> matrix{mat(0, 0) * scalar,
//                                                mat(0, 1) * scalar,
//                                                mat(0, 2) * scalar,
//                                                mat(1, 0) * scalar,
//                                                mat(1, 1) * scalar,
//                                                mat(1, 2) * scalar,
//                                                mat(2, 0) * scalar,
//                                                mat(2, 1) * scalar,
//                                                mat(2, 2) * scalar};
//   return matrix;
// }

// template <typename T1,
//           typename T2,
//           typename MetaPhysicL::boostcopy::enable_if_c<MetaPhysicL::ScalarTraits<T2>::value,
//                                                        int>::type = 0>
// auto operator*(const T2 & scalar, const TensorValue<T1> & mat)
//     -> TensorValue<decltype(mat(0, 0) * scalar)>
// {
//   return {mat(0, 0) * scalar,
//           mat(0, 1) * scalar,
//           mat(0, 2) * scalar,
//           mat(1, 0) * scalar,
//           mat(1, 1) * scalar,
//           mat(1, 2) * scalar,
//           mat(2, 0) * scalar,
//           mat(2, 1) * scalar,
//           mat(2, 2) * scalar};
// }

// template <typename T1,
//           typename T2,
//           typename MetaPhysicL::boostcopy::enable_if_c<MetaPhysicL::ScalarTraits<T2>::value,
//                                                        int>::type = 0>
// auto operator*(const TensorValue<T1> & mat, const T2 & scalar)
//     -> TensorValue<decltype(mat(0, 0) * scalar)>
// {

//   return {mat(0, 0) * scalar,
//           mat(0, 1) * scalar,
//           mat(0, 2) * scalar,
//           mat(1, 0) * scalar,
//           mat(1, 1) * scalar,
//           mat(1, 2) * scalar,
//           mat(2, 0) * scalar,
//           mat(2, 1) * scalar,
//           mat(2, 2) * scalar};
// }

// template <typename T1,
//           typename T2,
//           typename MetaPhysicL::boostcopy::enable_if_c<MetaPhysicL::ScalarTraits<T2>::value,
//                                                        int>::type = 0>
// auto operator*(const T2 & scalar, const TensorValue<T1> & mat)
//     -> TensorValue<decltype(mat(0, 0) * scalar)>
// {
//   return {mat(0, 0) * scalar,
//           mat(0, 1) * scalar,
//           mat(0, 2) * scalar,
//           mat(1, 0) * scalar,
//           mat(1, 1) * scalar,
//           mat(1, 2) * scalar,
//           mat(2, 0) * scalar,
//           mat(2, 1) * scalar,
//           mat(2, 2) * scalar};
// }

typedef DualNumber<VectorValue<Real>, NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<Real>>>
    ADRealVectorValue;
typedef ADRealVectorValue ADRealGradient;

typedef DualNumber<TensorValue<Real>, NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<Real>>>
    ADRealTensorValue;
typedef ADRealTensorValue ADRealTensor;

typedef MooseArray<ADReal> ADVariableValue;
typedef MooseArray<ADRealGradient> ADVariableGradient;
typedef MooseArray<ADRealTensor> ADVariableSecond;

#endif
