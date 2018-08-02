#ifndef MOOSEADTYPES_H
#define MOOSEADTYPES_H

#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/compare_types.h"

#include "metaphysicl/dualnumber.h"
#include "metaphysicl/numberarray.h"

template <typename>
class MooseArray;

// The 100 here is for how many DoFs there are per element.
#define AD_MAX_DOFS_PER_ELEM 100
using MetaPhysicL::DualNumber;
using MetaPhysicL::NumberArray;
using MetaPhysicL::TensorTraits;
using MetaPhysicL::VectorTraits;
using libMesh::Real;
using libMesh::ScalarTraits;
using libMesh::TensorValue;
using libMesh::TypeTensor;
using libMesh::TypeVector;
using libMesh::VectorValue;
using libMesh::boostcopy::enable_if_c;

/**
 * DualNumber naming
 */
template <typename T>
using ScalarDN = DualNumber<T, NumberArray<AD_MAX_DOFS_PER_ELEM, T>>;
template <typename T, template <class> class W>
using TemplateDN = DualNumber<W<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, W<T>>>;
template <typename T>
using VectorDN = DualNumber<VectorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>>>;
template <typename T>
using TensorDN = DualNumber<TensorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>>>;
template <typename T>
using VectorValueDN = DualNumber<VectorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>>>;
template <typename T>
using TensorValueDN = DualNumber<TensorValue<T>, NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>>>;

typedef ScalarDN<Real> ADReal;

// /*
//  * Vector-scalar multiplication
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TemplateDN<T, W> & vec, const T & scalar)
// {
//   VectorValue<T> value = vec.value() * scalar;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec.derivatives()[i] * scalar;
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const T & scalar, const TemplateDN<T, W> & vec)
// {
//   VectorValue<T> value = vec.value() * scalar;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec.derivatives()[i] * scalar;
//   return {value, derivatives};
// }

// template <typename T, typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TypeVector<T> & vec, const ScalarDN<T> & scalar)
// {
//   VectorValue<T> value = vec * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec * scalar.derivatives()[i];
//   return {value, derivatives};
// }

// template <typename T, typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const ScalarDN<T> & scalar, const TypeVector<T> & vec)
// {
//   VectorValue<T> value = vec * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec * scalar.derivatives()[i];
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TemplateDN<T, W> & vec, const ScalarDN<T> & scalar)
// {
//   VectorValue<T> value = vec.value() * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec.value() * scalar.derivatives()[i] + vec.derivatives()[i] *
//     scalar.value();
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const ScalarDN<T> & scalar, const TemplateDN<T, W> & vec)
// {
//   VectorValue<T> value = vec.value() * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = vec.value() * scalar.derivatives()[i] + vec.derivatives()[i] *
//     scalar.value();
//   return {value, derivatives};
// }

// /*
//  * Tensor-scalar multiplication
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const TemplateDN<T, W> & tensor, const T & scalar)
// {
//   TensorValue<T> value = tensor.value() * scalar;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor.derivatives()[i] * scalar;
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const T & scalar, const TemplateDN<T, W> & tensor)
// {
//   TensorValue<T> value = tensor.value() * scalar;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor.derivatives()[i] * scalar;
//   return {value, derivatives};
// }

// template <typename T, typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const TensorValue<T> & tensor, const ScalarDN<T> & scalar)
// {
//   TensorValue<T> value = tensor * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor * scalar.derivatives()[i];
//   return {value, derivatives};
// }

// template <typename T, typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const ScalarDN<T> & scalar, const TensorValue<T> & tensor)
// {
//   TensorValue<T> value = tensor * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor * scalar.derivatives()[i];
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const TemplateDN<T, W> & tensor, const ScalarDN<T> & scalar)
// {
//   TensorValue<T> value = tensor.value() * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] =
//         tensor.value() * scalar.derivatives()[i] + tensor.derivatives()[i] * scalar.value();
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T> operator*(const ScalarDN<T> & scalar, const TemplateDN<T, W> & tensor)
// {
//   TensorValue<T> value = tensor.value() * scalar.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] =
//         tensor.value() * scalar.derivatives()[i] + tensor.derivatives()[i] * scalar.value();
//   return {value, derivatives};
// }

// /*
//  * Vector-vector multiplication
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// ScalarDN<T> operator*(const TemplateDN<T, W> & dn, const TypeVector<T> & vec2)
// {
//   const VectorValue<T> & vec1 = dn.value();
//   T value = vec1 * vec2;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn.derivatives()[i] * vec2;
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// ScalarDN<T> operator*(const TypeVector<T> & vec2, const TemplateDN<T, W> & dn)
// {
//   const VectorValue<T> & vec1 = dn.value();
//   T value = vec1 * vec2;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn.derivatives()[i] * vec2;
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// ScalarDN<T> operator*(const TemplateDN<T, W> & dn1, const TemplateDN<T, W> & dn2)
// {
//   T value = dn1.value() * dn2.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, T> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn1.derivatives()[i] * dn2.value() + dn1.value() * dn2.derivatives()[i];
//   return {value, derivatives};
// }

// /*
//  * Vector-vector addition
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator+(const TemplateDN<T, W> & dn, const TypeVector<T> & vec2)
// {
//   return {dn.value() + vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator+(const TypeVector<T> & vec2, const TemplateDN<T, W> & dn)
// {
//   return {dn.value() + vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator+(const TemplateDN<T, W> & dn1, const TemplateDN<T, W> & dn2)
// {
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn1.derivatives()[i] + dn2.derivatives()[i];
//   return {dn1.value() + dn2.value(), derivatives};
// }

// /*
//  * Vector-vector subtraction
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator-(const TemplateDN<T, W> & dn, const TypeVector<T> & vec2)
// {
//   return {dn.value() - vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator-(const TypeVector<T> & vec2, const TemplateDN<T, W> & dn)
// {
//   return {vec2 - dn.value(), -dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T>
// operator-(const TemplateDN<T, W> & dn1, const TemplateDN<T, W> & dn2)
// {
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn1.derivatives()[i] - dn2.derivatives()[i];
//   return {dn1.value() - dn2.value(), derivatives};
// }

// /*
//  * Tensor-tensor addition
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator+(const TemplateDN<T, W> & dn, const TypeTensor<T> & vec2)
// {
//   return {dn.value() + vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator+(const TypeTensor<T> & vec2, const TemplateDN<T, W> & dn)
// {
//   return {dn.value() + vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator+(const TemplateDN<T, W> & dn1, const TemplateDN<T, W> & dn2)
// {
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn1.derivatives()[i] + dn2.derivatives()[i];
//   return {dn1.value() + dn2.value(), derivatives};
// }

// /*
//  * Tensor-tensor subtraction
//  */
// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator-(const TemplateDN<T, W> & dn, const TypeTensor<T> & vec2)
// {
//   return {dn.value() - vec2, dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator-(const TypeTensor<T> & vec2, const TemplateDN<T, W> & dn)
// {
//   return {vec2 - dn.value(), -dn.derivatives()};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// TensorDN<T>
// operator-(const TemplateDN<T, W> & dn1, const TemplateDN<T, W> & dn2)
// {
//   NumberArray<AD_MAX_DOFS_PER_ELEM, TensorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = dn1.derivatives()[i] - dn2.derivatives()[i];
//   return {dn1.value() - dn2.value(), derivatives};
// }

// /*
//  * Tensor-vector multiplication
//  */
// template <template <typename> class W,
//           template <typename> class W2,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<VectorTraits<W2, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TemplateDN<T, W> & tensor, const TemplateDN<T, W2> & vector)
// {
//   VectorValue<T> value = tensor.value() * vector.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] =
//         tensor.value() * vector.derivatives()[i] + tensor.derivatives()[i] * vector.value();
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<TensorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TemplateDN<T, W> & tensor, const TypeVector<T> & vector)
// {
//   VectorValue<T> value = tensor.value() * vector;
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor.derivatives()[i] * vector;
//   return {value, derivatives};
// }

// template <template <typename> class W,
//           typename T,
//           typename enable_if_c<VectorTraits<W, T>::value, int>::type = 0,
//           typename enable_if_c<ScalarTraits<T>::value, int>::type = 0>
// VectorDN<T> operator*(const TypeTensor<T> & tensor, const TemplateDN<T, W> & vector)
// {
//   VectorValue<T> value = tensor * vector.value();
//   NumberArray<AD_MAX_DOFS_PER_ELEM, VectorValue<T>> derivatives;
//   for (decltype(AD_MAX_DOFS_PER_ELEM) i = 0; i < AD_MAX_DOFS_PER_ELEM; ++i)
//     derivatives[i] = tensor * vector.derivatives()[i];
//   return {value, derivatives};
// }

/*
 * Some helpful typedefs
 */
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
