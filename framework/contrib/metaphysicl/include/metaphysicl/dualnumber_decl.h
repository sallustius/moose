//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MetaPhysicL - A metaprogramming library for physics calculations
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------

#ifndef METAPHYSICL_DUALNUMBER_DECL_H
#define METAPHYSICL_DUALNUMBER_DECL_H

#include <ostream>
#include <limits>
#include <utility>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/dualderivatives.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/testable.h"
#include "metaphysicl/numberarray.h"

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "MooseError.h"

namespace MetaPhysicL
{
template <typename T, size_t N, typename Enable = void>
class DualNumber;

template <typename T, typename D = T>
class DualNumberBase : public safe_bool<DualNumberBase<T, D>>
{
public:
  typedef T value_type;

  typedef D derivatives_type;

  DualNumberBase();

  template <typename T2>
  DualNumberBase(const T2 & val);

  template <typename T2, typename D2>
  DualNumberBase(const T2 & val, const D2 & deriv);

#if __cplusplus >= 201103L
  // Move constructors are useful when all your data is on the heap
  DualNumberBase(DualNumberBase<T, D> && src) = default;

  // Move assignment avoids heap operations too
  DualNumberBase & operator=(DualNumberBase<T, D> && src) = default;

  // Standard copy operations get implicitly deleted upon move
  // constructor definition, so we redefine them.
  DualNumberBase(const DualNumberBase<T, D> & src) = default;

  DualNumberBase & operator=(const DualNumberBase<T, D> & src) = default;
#endif

  template <typename T2>
  DualNumberBase & operator=(const T2 & value)
  {
    _val = value;
    _deriv.zero();
    return *this;
  }

  T & value();

  const T & value() const;

  D & derivatives();

  const D & derivatives() const;

  bool boolean_test() const;

  DualNumberBase<T, D> operator-() const;

  DualNumberBase<T, D> operator!() const;

  template <typename T2, typename D2>
  auto operator+=(const DualNumberBase<T2, D2> & a) -> decltype(*this) &;

  template <typename T2, size_t N>
  auto operator+=(const DualNumber<T2, N> & a) -> decltype(*this) &;

  template <typename T2>
  auto operator+=(const T2 & a) -> decltype(*this) &;

  template <typename T2, typename D2>
  auto operator-=(const DualNumberBase<T2, D2> & a) -> decltype(*this) &;

  template <typename T2, size_t N>
  auto operator-=(const DualNumber<T2, N> & a) -> decltype(*this) &;

  template <typename T2>
  auto operator-=(const T2 & a) -> decltype(*this) &;

  template <typename T2, typename D2>
  auto operator*=(const DualNumberBase<T2, D2> & a) -> decltype(*this) &;

  template <typename T2, size_t N>
  auto operator*=(const DualNumber<T2, N> & a) -> decltype(*this) &;

  template <typename T2>
  auto operator*=(const T2 & a) -> decltype(*this) &;

  template <typename T2, typename D2>
  auto operator/=(const DualNumberBase<T2, D2> & a) -> decltype(*this) &;

  template <typename T2, size_t N>
  auto operator/=(const DualNumber<T2, N> & a) -> decltype(*this) &;

  template <typename T2>
  auto operator/=(const T2 & a) -> decltype(*this) &;

  operator T &()
  {
    mooseDoOnce(mooseWarning("This message will only appear once: You are operating on a "
                             "DualNumber value without operating on its derivatives. "
                             "This may have a deleterious impact on your Jacobian calculations!"));
    return _val;
  }
  operator const T &() const { return _val; }

protected:
  T _val;
  D _deriv;
};

template <typename T, std::size_t N>
struct DualNumberSurrogate;

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumber<R, N> return_dn(const R & (TBase::*fn)(PtrArgs...) const,
                           const DualNumber<T, N> & calling_dn,
                           ParamArgs &&... args);

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumber<R, N> return_dn(R (TBase::*fn)(PtrArgs...) const,
                           const DualNumber<T, N> & calling_dn,
                           ParamArgs &&... args);

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumberSurrogate<R, N> &
return_dns(R & (TBase::*fn)(PtrArgs...), DualNumber<T, N> & calling_dn, ParamArgs &&... args);

template <typename T, typename TBase, size_t N, class... PtrArgs, class... ParamArgs>
void const_void_helper(void (TBase::*fn)(PtrArgs...) const,
                       const DualNumber<T, N> & calling_dn,
                       ParamArgs &&... args);

template <typename T, typename TBase, size_t N, class... PtrArgs, class... ParamArgs>
void
void_helper(void (TBase::*fn)(PtrArgs...), DualNumber<T, N> & calling_dn, ParamArgs &&... args);

// Base template
template <typename T, std::size_t N, typename Enable>
class DualNumber : public DualNumberBase<T, NumberArray<N, T>>
{
public:
  using DualNumberBase<T, NumberArray<N, T>>::DualNumberBase;

  DualNumber<T, N> operator-() const { return DualNumber<T, N>(-this->_val, -this->_deriv); }
  DualNumber<T, N> operator!() const { return DualNumber<T, N>(!this->_val, !this->_deriv); }
};

#define ConstReturnDecl(methodName)                                                                \
  template <class... Args>                                                                         \
  auto methodName(Args &&... args)                                                                 \
      const->DualNumber<typename std::remove_const<typename std::remove_reference<decltype(        \
                            this->value().methodName(std::forward<Args>(args)...))>::type>::type,  \
                        N>

#define ConstReturnDef(methodName, condition)                                                      \
  template <typename T, std::size_t N>                                                             \
  template <class... Args>                                                                         \
  auto DualNumber<T, N, typename std::enable_if<condition>::type>::methodName(Args &&... args)     \
      const->DualNumber<typename std::remove_const<typename std::remove_reference<decltype(        \
                            this->value().methodName(std::forward<Args>(args)...))>::type>::type,  \
                        N>                                                                         \
  {                                                                                                \
    return return_dn(&T::methodName, *this, std::forward<Args>(args)...);                          \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

#define NonConstReturnDecl(methodName)                                                             \
  template <class... Args>                                                                         \
  auto methodName(Args &&... args)                                                                 \
      ->DualNumberSurrogate<typename std::remove_reference<decltype(                               \
                                this->value().methodName(std::forward<Args>(args)...))>::type,     \
                            N> &

#define NonConstReturnDef(methodName, condition)                                                   \
  template <typename T, std::size_t N>                                                             \
  template <class... Args>                                                                         \
  auto DualNumber<T, N, typename std::enable_if<condition>::type>::methodName(Args &&... args)     \
      ->DualNumberSurrogate<typename std::remove_reference<decltype(                               \
                                this->value().methodName(std::forward<Args>(args)...))>::type,     \
                            N> &                                                                   \
  {                                                                                                \
    return return_dns(&T::methodName, *this, std::forward<Args>(args)...);                         \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

#define ConstVoidDecl(methodName)                                                                  \
  template <class... Args>                                                                         \
  void methodName(Args &&... args) const

#define ConstVoidDef(methodName, condition)                                                        \
  template <typename T, std::size_t N>                                                             \
  template <class... Args>                                                                         \
  void DualNumber<T, N, typename std::enable_if<condition>::type>::methodName(Args &&... args)     \
      const                                                                                        \
  {                                                                                                \
    const_void_helper(&T::methodName, *this, std::forward<Args>(args)...);                         \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

#define NonConstVoidDecl(methodName)                                                               \
  template <class... Args>                                                                         \
  void methodName(Args &&... args)

#define NonConstVoidDef(methodName, condition)                                                     \
  template <typename T, std::size_t N>                                                             \
  template <class... Args>                                                                         \
  void DualNumber<T, N, typename std::enable_if<condition>::type>::methodName(Args &&... args)     \
                                                                                                   \
  {                                                                                                \
    void_helper(&T::methodName, *this, std::forward<Args>(args)...);                               \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

// Vector
template <typename T, std::size_t N>
class DualNumber<T,
                 N,
                 typename std::enable_if<is_same_template<T, TypeVector>::value ||
                                         is_same_template<T, VectorValue>::value>::type>
  : public DualNumberBase<T, NumberArray<N, T>>
{
public:
  using DualNumberBase<T, NumberArray<N, T>>::DualNumberBase;

  DualNumber<T, N> operator-() const { return DualNumber<T, N>(-this->_val, -this->_deriv); }
  DualNumber<T, N> operator!() const { return DualNumber<T, N>(!this->_val, !this->_deriv); }

  ConstReturnDecl(operator());
  NonConstReturnDecl(operator());
  NonConstVoidDecl(zero);

  std::map<typename T::index_type, DualNumberSurrogate<typename T::value_type, N>>
      dual_number_surrogates;
};

#ifndef COMMA
#define COMMA ,
#endif

ConstReturnDef(operator(),
               is_same_template<T COMMA TypeVector>::value ||
                   is_same_template<T COMMA VectorValue>::value);
NonConstReturnDef(operator(),
                  is_same_template<T COMMA TypeVector>::value ||
                      is_same_template<T COMMA VectorValue>::value);
NonConstVoidDef(zero,
                is_same_template<T COMMA TypeVector>::value ||
                    is_same_template<T COMMA VectorValue>::value);

// Tensor
template <typename T, std::size_t N>
class DualNumber<T,
                 N,
                 typename std::enable_if<is_same_template<T, TypeTensor>::value ||
                                         is_same_template<T, TensorValue>::value>::type>
  : public DualNumberBase<T, NumberArray<N, T>>
{
public:
  using DualNumberBase<T, NumberArray<N, T>>::DualNumberBase;

  DualNumber<T, N> operator-() const { return DualNumber<T, N>(-this->_val, -this->_deriv); }
  DualNumber<T, N> operator!() const { return DualNumber<T, N>(!this->_val, !this->_deriv); }

  ConstReturnDecl(operator());
  NonConstReturnDecl(operator());
  NonConstVoidDecl(zero);
  ConstReturnDecl(row);
  ConstReturnDecl(tr);
  ConstReturnDecl(transpose);
  NonConstVoidDecl(rotate);

  std::map<typename T::index_type, DualNumberSurrogate<typename T::value_type, N>>
      dual_number_surrogates;
};

ConstReturnDef(operator(),
               is_same_template<T COMMA TypeTensor>::value ||
                   is_same_template<T COMMA TensorValue>::value);
NonConstReturnDef(operator(),
                  is_same_template<T COMMA TypeTensor>::value ||
                      is_same_template<T COMMA TensorValue>::value);
NonConstVoidDef(zero,
                is_same_template<T COMMA TypeTensor>::value ||
                    is_same_template<T COMMA TensorValue>::value);
ConstReturnDef(row,
               is_same_template<T COMMA TypeTensor>::value ||
                   is_same_template<T COMMA TensorValue>::value);
ConstReturnDef(tr,
               is_same_template<T COMMA TypeTensor>::value ||
                   is_same_template<T COMMA TensorValue>::value);
ConstReturnDef(transpose,
               is_same_template<T COMMA TypeTensor>::value ||
                   is_same_template<T COMMA TensorValue>::value);
NonConstVoidDef(rotate,
                is_same_template<T COMMA TypeTensor>::value ||
                    is_same_template<T COMMA TensorValue>::value);

// RankTwo
template <typename T, std::size_t N>
class DualNumber<T, N, typename std::enable_if<std::is_same<T, RankTwoTensor>::value>::type>
  : public DualNumberBase<T, NumberArray<N, T>>
{
public:
  using DualNumberBase<T, NumberArray<N, T>>::DualNumberBase;

  DualNumber<T, N> operator-() const { return DualNumber<T, N>(-this->_val, -this->_deriv); }
  DualNumber<T, N> operator!() const { return DualNumber<T, N>(!this->_val, !this->_deriv); }

  template <typename T2>
  DualNumber(const DualNumber<TypeVector<T2>, N> & row1,
             const DualNumber<TypeVector<T2>, N> & row2,
             const DualNumber<TypeVector<T2>, N> & row3)
    : DualNumberBase<T, NumberArray<N, T>>()
  {
    this->value() = RankTwoTensor(row1.value(), row2.value(), row3.value());
    for (decltype(N) i = 0; i < N; ++i)
      this->derivatives()[i] =
          RankTwoTensor(row1.derivatives()[i], row2.derivatives()[i], row3.derivatives()[i]);
  }
  template <typename T2>
  DualNumber(const DualNumber<VectorValue<T2>, N> & row1,
             const DualNumber<VectorValue<T2>, N> & row2,
             const DualNumber<VectorValue<T2>, N> & row3)
    : DualNumberBase<T, NumberArray<N, T>>()
  {
    this->value() = RankTwoTensor(row1.value(), row2.value(), row3.value());
    for (decltype(N) i = 0; i < N; ++i)
      this->derivatives()[i] =
          RankTwoTensor(row1.derivatives()[i], row2.derivatives()[i], row3.derivatives()[i]);
  }

  ConstReturnDecl(operator());
  NonConstReturnDecl(operator());
  NonConstVoidDecl(zero);
  ConstReturnDecl(row);
  ConstReturnDecl(tr);
  ConstReturnDecl(transpose);
  NonConstVoidDecl(rotate);
  ConstReturnDecl(positiveProjectionEigenDecomposition);
  ConstVoidDecl(symmetricEigenvaluesEigenvectors);
  ConstReturnDecl(deviatoric);

  std::map<typename T::index_type, DualNumberSurrogate<typename T::value_type, N>>
      dual_number_surrogates;
};

ConstReturnDef(operator(), std::is_same<T COMMA RankTwoTensor>::value);
NonConstReturnDef(operator(), std::is_same<T COMMA RankTwoTensor>::value);
NonConstVoidDef(zero, std::is_same<T COMMA RankTwoTensor>::value);
ConstReturnDef(row, std::is_same<T COMMA RankTwoTensor>::value);
ConstReturnDef(tr, std::is_same<T COMMA RankTwoTensor>::value);
ConstReturnDef(transpose, std::is_same<T COMMA RankTwoTensor>::value);
NonConstVoidDef(rotate, std::is_same<T COMMA RankTwoTensor>::value);
ConstReturnDef(positiveProjectionEigenDecomposition, std::is_same<T COMMA RankTwoTensor>::value);
ConstVoidDef(symmetricEigenvaluesEigenvectors, std::is_same<T COMMA RankTwoTensor>::value);
ConstReturnDef(deviatoric, std::is_same<T COMMA RankTwoTensor>::value);

// RankFour
template <typename T, std::size_t N>
class DualNumber<T, N, typename std::enable_if<std::is_same<T, RankFourTensor>::value>::type>
  : public DualNumberBase<T, NumberArray<N, T>>
{
public:
  using DualNumberBase<T, NumberArray<N, T>>::DualNumberBase;
  using DualNumberBase<T, NumberArray<N, T>>::operator=;

  DualNumber<T, N> operator-() const { return DualNumber<T, N>(-this->_val, -this->_deriv); }
  DualNumber<T, N> operator!() const { return DualNumber<T, N>(!this->_val, !this->_deriv); }

  ConstReturnDecl(operator());
  NonConstReturnDecl(operator());
  NonConstVoidDecl(zero);
  NonConstVoidDecl(rotate);
  ConstReturnDecl(invSymm);

  std::map<typename T::index_type, DualNumberSurrogate<Real, N>> dual_number_surrogates;
};

ConstReturnDef(operator(), std::is_same<T COMMA RankFourTensor>::value);
NonConstReturnDef(operator(), std::is_same<T COMMA RankFourTensor>::value);
NonConstVoidDef(zero, std::is_same<T COMMA RankFourTensor>::value);
NonConstVoidDef(rotate, std::is_same<T COMMA RankFourTensor>::value);
ConstReturnDef(invSymm, std::is_same<T COMMA RankFourTensor>::value);

template <typename T, std::size_t N>
struct DualNumberSurrogate
{
  DualNumberSurrogate(DualNumber<T, N> & dn);
  DualNumberSurrogate(DualNumber<T, N> && dn);
  DualNumberSurrogate(const T & n);

  template <typename T2, class... Args>
  DualNumberSurrogate(DualNumber<T2, N> & dn, Args &&... args)
    : value(dn.value()(std::forward<Args>(args)...))
  {
    for (decltype(N) di = 0; di < N; ++di)
      derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
  }

  template <typename T2, class... Args>
  DualNumberSurrogate(DualNumber<T2, N> && dn, Args &&... args)
    : value(dn.value()(std::forward<Args>(args)...))
  {
    for (decltype(N) di = 0; di < N; ++di)
      derivatives[di] = &dn.derivatives()[di](std::forward<Args>(args)...);
  }

  DualNumberSurrogate(DualNumberSurrogate<T, N> & dns)
    : value(dns.value), derivatives(dns.derivatives)
  {
  }

  DualNumberSurrogate(const DualNumberSurrogate<T, N> & dns)
    : value(dns.value), derivatives(dns.derivatives)
  {
  }

  DualNumberSurrogate(DualNumberSurrogate<T, N> && dns)
    : value(dns.value), derivatives(dns.derivatives)
  {
  }

  DualNumberSurrogate<T, N> & operator=(DualNumberSurrogate<T, N> & dns)
  {
    value = dns.value;
    for (decltype(N) i = 0; i < N; ++i)
    {
      if (dns.derivatives[i])
        *derivatives[i] = *dns.derivatives[i];
      else
        *derivatives[i] = 0;
    }
    return *this;
  }
  DualNumberSurrogate<T, N> & operator=(const DualNumberSurrogate<T, N> & dns)
  {
    value = dns.value;
    for (decltype(N) i = 0; i < N; ++i)
    {
      if (dns.derivatives[i])
        *derivatives[i] = *dns.derivatives[i];
      else
        *derivatives[i] = 0;
    }
    return *this;
  }
  DualNumberSurrogate<T, N> & operator=(DualNumberSurrogate<T, N> && dns)
  {
    value = dns.value;
    for (decltype(N) i = 0; i < N; ++i)
    {
      if (dns.derivatives[i])
        *derivatives[i] = *dns.derivatives[i];
      else
        *derivatives[i] = 0;
    }
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator=(const T2 & in_value)
  {
    value = in_value;
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] = 0;
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator=(const DualNumber<T2, N> & in_dn)
  {
    value = in_dn.value();
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] = in_dn.derivatives()[i];
    return *this;
  }

  DualNumberSurrogate<T, N> & operator+=(const DualNumberSurrogate<T, N> & dns)
  {
    value += dns.value;
    for (decltype(N) i = 0; i < N; ++i)
      if (dns.derivatives[i])
        *derivatives[i] += *dns.derivatives[i];
    return *this;
  }
  DualNumberSurrogate<T, N> & operator-=(const DualNumberSurrogate<T, N> & dns)
  {
    value -= dns.value;
    for (decltype(N) i = 0; i < N; ++i)
      if (dns.derivatives[i])
        *derivatives[i] -= *dns.derivatives[i];
    return *this;
  }
  DualNumberSurrogate<T, N> & operator*=(const DualNumberSurrogate<T, N> & dns)
  {
    value *= dns.value;
    for (decltype(N) i = 0; i < N; ++i)
      if (dns.derivatives[i])
        *derivatives[i] *= *dns.derivatives[i];
    return *this;
  }
  DualNumberSurrogate<T, N> & operator/=(const DualNumberSurrogate<T, N> & dns)
  {
    value /= dns.value;
    for (decltype(N) i = 0; i < N; ++i)
      if (dns.derivatives[i])
        *derivatives[i] /= *dns.derivatives[i];
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator+=(const T2 & in_value)
  {
    value += in_value;
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator-=(const T2 & in_value)
  {
    value -= in_value;
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator*=(const T2 & in_value)
  {
    value *= in_value;
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator/=(const T2 & in_value)
  {
    value /= in_value;
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator+=(const DualNumber<T, N> & in_dn)
  {
    value += in_dn.value();
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] += in_dn.derivatives()[i];
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator-=(const DualNumber<T, N> & in_dn)
  {
    value -= in_dn.value();
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] -= in_dn.derivatives()[i];
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator*=(const DualNumber<T, N> & in_dn)
  {
    value *= in_dn.value();
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] *= in_dn.derivatives()[i];
    return *this;
  }
  template <typename T2>
  DualNumberSurrogate<T, N> & operator/=(const DualNumber<T, N> & in_dn)
  {
    value /= in_dn.value();
    for (decltype(N) i = 0; i < N; ++i)
      *derivatives[i] /= in_dn.derivatives()[i];
    return *this;
  }

  operator T() const { return value; }

  T & value;
  NumberArray<N, T *> derivatives;
};

#define DNS_compares(comparator)                                                                   \
  template <typename T, size_t N>                                                                  \
  bool operator comparator(const DualNumberSurrogate<T, N> & dns, const T & solo)                  \
  {                                                                                                \
    return dns.value comparator solo;                                                              \
  }                                                                                                \
  template <typename T, size_t N>                                                                  \
  bool operator comparator(const T & solo, const DualNumberSurrogate<T, N> & dns)                  \
  {                                                                                                \
    return solo comparator dns.value;                                                              \
  }                                                                                                \
  template <typename T, size_t N>                                                                  \
  bool operator comparator(const DualNumberSurrogate<T, N> & dns1,                                 \
                           const DualNumberSurrogate<T, N> & dns2)                                 \
  {                                                                                                \
    return dns1.value comparator dns2.value;                                                       \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

DNS_compares(>);
DNS_compares(<);
DNS_compares(==);
DNS_compares(!=);
DNS_compares(>=);
DNS_compares(<=);

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumber<R, N>
return_dn(const R & (TBase::*fn)(PtrArgs...) const,
          const DualNumber<T, N> & calling_dn,
          ParamArgs &&... args)
{
  NumberArray<N, R> deriv;
  for (decltype(N) di = 0; di < N; ++di)
    deriv[di] = (calling_dn.derivatives()[di].*fn)(std::forward<ParamArgs>(args)...);
  return {(calling_dn.value().*fn)(std::forward<ParamArgs>(args)...), deriv};
}

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumber<R, N>
return_dn(R (TBase::*fn)(PtrArgs...) const,
          const DualNumber<T, N> & calling_dn,
          ParamArgs &&... args)
{
  NumberArray<N, R> deriv;
  for (decltype(N) di = 0; di < N; ++di)
    deriv[di] = (calling_dn.derivatives()[di].*fn)(std::forward<ParamArgs>(args)...);
  return {(calling_dn.value().*fn)(std::forward<ParamArgs>(args)...), deriv};
}

template <typename T, typename TBase, size_t N, typename R, class... PtrArgs, class... ParamArgs>
DualNumberSurrogate<R, N> &
return_dns(R & (TBase::*)(PtrArgs...), DualNumber<T, N> & calling_dn, ParamArgs &&... args)
{
  std::tuple<PtrArgs...> key(std::forward<ParamArgs>(args)...);
  typename std::map<decltype(key), DualNumberSurrogate<R, N>>::iterator it =
      calling_dn.dual_number_surrogates.find(key);
  if (it == calling_dn.dual_number_surrogates.end())
  {
    DualNumberSurrogate<R, N> dns{calling_dn, std::forward<ParamArgs>(args)...};
    calling_dn.dual_number_surrogates.emplace(key, std::move(dns));
  }
  return calling_dn.dual_number_surrogates.at(key);
}

template <typename T, typename TBase, size_t N, class... PtrArgs, class... ParamArgs>
void
const_void_helper(void (TBase::*fn)(PtrArgs...) const,
                  const DualNumber<T, N> & calling_dn,
                  ParamArgs &&... args)
{
  (calling_dn.value().*fn)(std::forward<ParamArgs>(args)...);
  for (decltype(N) i = 0; i < N; ++i)
    (calling_dn.derivatives()[i].*fn)(std::forward<ParamArgs>(args)...);
}

template <typename T, typename TBase, size_t N, class... PtrArgs, class... ParamArgs>
void
void_helper(void (TBase::*fn)(PtrArgs...), DualNumber<T, N> & calling_dn, ParamArgs &&... args)
{
  (calling_dn.value().*fn)(std::forward<ParamArgs>(args)...);
  for (decltype(N) i = 0; i < N; ++i)
    (calling_dn.derivatives()[i].*fn)(std::forward<ParamArgs>(args)...);
}

// Helper class to handle partial specialization for DualNumberBase
// constructors

template <typename T, typename D>
struct DualNumberConstructor
{
  static T value(const DualNumberBase<T, D> & v) { return v.value(); }

  template <typename T2>
  static T value(const T2 & v)
  {
    return v;
  }

  template <typename T2, typename D2>
  static T value(const T2 & v, const D2 &)
  {
    return v;
  }

  template <typename T2, typename D2>
  static T value(const DualNumberBase<T2, D2> & v)
  {
    return DualNumberConstructor<T, D>::value(v.value());
  }

  template <typename T2, size_t N>
  static T value(const DualNumber<T2, N> & v)
  {
    return DualNumberConstructor<T, NumberArray<N, T>>::value(v.value());
  }

  template <typename T2>
  static D deriv(const T2 &)
  {
    return {};
  }

  template <typename T2, typename D2>
  static D deriv(const DualNumberBase<T2, D2> & v)
  {
    return v.derivatives();
  }

  template <typename T2, size_t N>
  static D deriv(const DualNumber<T2, N> & v)
  {
    return v.derivatives();
  }

  template <typename T2, typename D2>
  static D deriv(const T2 &, const D2 & d)
  {
    return d;
  }
};

template <typename T, typename D, typename DD>
struct DualNumberConstructor<DualNumberBase<T, D>, DD>
{
  template <typename T2, typename D2, typename D3>
  static DualNumberBase<T, D> value(const DualNumberBase<DualNumberBase<T2, D2>, D3> & v)
  {
    return v.value();
  }

  template <typename T2>
  static DualNumberBase<T, D> value(const T2 & v)
  {
    return v;
  }

  template <typename T2, typename D2>
  static DualNumberBase<T, D> value(const T2 & v, const D2 & d)
  {
    return DualNumberBase<T, D>(v, d);
  }

  template <typename D2>
  static DualNumberBase<T, D> value(const DualNumberBase<T, D> & v, const D2 &)
  {
    return v;
  }

  template <typename T2>
  static DD deriv(const T2 &)
  {
    return 0;
  }

  template <typename T2, typename D2>
  static DD deriv(const DualNumberBase<T2, D2> & v)
  {
    return v.derivatives();
  }

  template <typename T2, typename D2>
  static DD deriv(const T2 &, const D2 & d)
  {
    return d;
  }
};

// FIXME: these operators currently do automatic type promotion when
// encountering DualNumbers of differing levels of recursion and
// differentiability.  But what we really want is automatic type
// *demotion*, to avoid pretending we have accurate derivatives which
// we don't have.  If we could do that right then some potential
// subtle run-time user errors would turn into compile-time user
// errors.

#define DualNumber_decl_op(opname, functorname)                                                    \
  template <typename T, size_t N, typename T2>                                                     \
  inline auto operator opname(const DualNumber<T, N> & a, const DualNumber<T2, N> & b)             \
      ->DualNumber<decltype(a.value() opname b.value()), N>;                                       \
                                                                                                   \
  template <typename T, typename T2, size_t N>                                                     \
  inline auto operator opname(const T & a, const DualNumber<T2, N> & b)                            \
      ->DualNumber<decltype(a opname b.value()), N>;                                               \
                                                                                                   \
  template <typename T, size_t N, typename T2>                                                     \
  inline auto operator opname(const DualNumber<T, N> & a, const T2 & b)                            \
      ->DualNumber<decltype(a.value() opname b), N>

DualNumber_decl_op(+, Plus);
DualNumber_decl_op(-, Minus);
DualNumber_decl_op(*, Multiplies);
DualNumber_decl_op(/, Divides);

#define DualNumber_decl_compare(opname)                                                            \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline bool operator opname(const DualNumberBase<T, D> & a, const DualNumberBase<T2, D2> & b);   \
                                                                                                   \
  template <typename T, typename T2, typename D2>                                                  \
  inline typename boostcopy::                                                                      \
      enable_if_class<typename CompareTypes<DualNumberBase<T2, D2>, T>::supertype, bool>::type     \
      operator opname(const T & a, const DualNumberBase<T2, D2> & b);                              \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename boostcopy::                                                                      \
      enable_if_class<typename CompareTypes<DualNumberBase<T, D>, T2>::supertype, bool>::type      \
      operator opname(const DualNumberBase<T, D> & a, const T2 & b)

DualNumber_decl_compare(>);
DualNumber_decl_compare(>=);
DualNumber_decl_compare(<);
DualNumber_decl_compare(<=);
DualNumber_decl_compare(==);
DualNumber_decl_compare(!=);
DualNumber_decl_compare(&&);
DualNumber_decl_compare(||);

template <typename T, typename D>
inline std::ostream & operator<<(std::ostream & output, const DualNumberBase<T, D> & a);

// ScalarTraits, RawType, CompareTypes specializations

template <typename T, typename D>
struct ScalarTraits<DualNumberBase<T, D>>
{
  static const bool value = ScalarTraits<T>::value;
};

template <typename T, typename D>
struct RawType<DualNumberBase<T, D>>
{
  typedef typename RawType<T>::value_type value_type;

  static value_type value(const DualNumberBase<T, D> & a) { return raw_value(a.value()); }
};

template <typename T, typename T2, typename D, bool reverseorder>
struct PlusType<DualNumberBase<T, D>,
                T2,
                reverseorder,
                typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricPlusType<T, T2, reverseorder>::supertype, D> supertype;
};

template <typename T, typename D, typename T2, typename D2, bool reverseorder>
struct PlusType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricPlusType<T, T2, reverseorder>::supertype,
                         typename SymmetricPlusType<D, D2, reverseorder>::supertype>
      supertype;
};

template <typename T, typename D, bool reverseorder>
struct PlusType<DualNumberBase<T, D>, DualNumberBase<T, D>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricPlusType<T, T>::supertype,
                         typename SymmetricPlusType<D, D>::supertype>
      supertype;
};

template <typename T, typename T2, typename D, bool reverseorder>
struct MinusType<DualNumberBase<T, D>,
                 T2,
                 reverseorder,
                 typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricMinusType<T, T2, reverseorder>::supertype, D> supertype;
};

template <typename T, typename D, typename T2, typename D2, bool reverseorder>
struct MinusType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricMinusType<T, T2, reverseorder>::supertype,
                         typename SymmetricMinusType<D, D2, reverseorder>::supertype>
      supertype;
};

template <typename T, typename D, bool reverseorder>
struct MinusType<DualNumberBase<T, D>, DualNumberBase<T, D>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricMinusType<T, T>::supertype,
                         typename SymmetricMinusType<D, D>::supertype>
      supertype;
};

template <typename T, typename T2, typename D, bool reverseorder>
struct MultipliesType<DualNumberBase<T, D>, T2, reverseorder>
// typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricMultipliesType<T, T2, reverseorder>::supertype,
                         typename SymmetricMultipliesType<D, T2, reverseorder>::supertype>
      supertype;
};

template <typename T, typename D, typename T2, typename D2, bool reverseorder>
struct MultipliesType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, reverseorder>
{
  typedef DualNumberBase<
      typename SymmetricMultipliesType<T, T2, reverseorder>::supertype,
      typename SymmetricPlusType<
          typename SymmetricMultipliesType<T, D2, reverseorder>::supertype,
          typename SymmetricMultipliesType<D, T2, reverseorder>::supertype>::supertype>
      supertype;
};

template <typename T, typename D, bool reverseorder>
struct MultipliesType<DualNumberBase<T, D>, DualNumberBase<T, D>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricMultipliesType<T, T, reverseorder>::supertype,
                         typename SymmetricMultipliesType<T, D, reverseorder>::supertype>
      supertype;
};

template <typename T, typename T2, typename D>
struct DividesType<DualNumberBase<T, D>,
                   T2,
                   false,
                   typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricDividesType<T, T2>::supertype,
                         typename SymmetricDividesType<D, T2>::supertype>
      supertype;
};

template <typename T, typename D, typename T2>
struct DividesType<DualNumberBase<T, D>,
                   T2,
                   true,
                   typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<
      typename SymmetricDividesType<T2, T>::supertype,
      typename SymmetricDividesType<typename SymmetricMultipliesType<T2, D>::supertype,
                                    T>::supertype>
      supertype;
};

template <typename T, typename D, typename T2, typename D2>
struct DividesType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, false>
{
  typedef DualNumberBase<
      typename SymmetricDividesType<T, T2>::supertype,
      typename SymmetricMinusType<
          typename SymmetricDividesType<T2, D>::supertype,
          typename SymmetricDividesType<typename SymmetricMultipliesType<T, D2>::supertype,
                                        T2>::supertype>::supertype>
      supertype;
};

template <typename T, typename D, typename T2, typename D2>
struct DividesType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, true>
{
  typedef typename DividesType<DualNumberBase<T2, D2>, DualNumberBase<T, D>, false>::supertype
      supertype;
};

template <typename T, typename D>
struct DividesType<DualNumberBase<T, D>, DualNumberBase<T, D>, false>
{
  typedef DualNumberBase<
      T,
      typename SymmetricMinusType<
          typename SymmetricDividesType<T, D>::supertype,
          typename SymmetricDividesType<typename SymmetricMultipliesType<T, D>::supertype,
                                        T>::supertype>::supertype>
      supertype;
};

template <typename T, typename D>
struct DividesType<DualNumberBase<T, D>, DualNumberBase<T, D>, true>
{
  typedef
      typename DividesType<DualNumberBase<T, D>, DualNumberBase<T, D>, false>::supertype supertype;
};

template <typename T, typename T2, typename D, bool reverseorder>
struct AndType<DualNumberBase<T, D>,
               T2,
               reverseorder,
               typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricAndType<T, T2, reverseorder>::supertype, bool> supertype;
};

template <typename T, typename D, typename T2, typename D2, bool reverseorder>
struct AndType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricAndType<T, T2, reverseorder>::supertype, bool> supertype;
};

template <typename T, typename D>
struct AndType<DualNumberBase<T, D>, DualNumberBase<T, D>>
{
  typedef DualNumberBase<typename SymmetricAndType<T, T>::supertype, bool> supertype;
};

template <typename T, typename T2, typename D, bool reverseorder>
struct OrType<DualNumberBase<T, D>,
              T2,
              reverseorder,
              typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<typename SymmetricOrType<T, T2, reverseorder>::supertype, bool> supertype;
};

template <typename T, typename D, typename T2, typename D2, bool reverseorder>
struct OrType<DualNumberBase<T, D>, DualNumberBase<T2, D2>, reverseorder>
{
  typedef DualNumberBase<typename SymmetricOrType<T, T2, reverseorder>::supertype, bool> supertype;
};

template <typename T, typename D>
struct OrType<DualNumberBase<T, D>, DualNumberBase<T, D>>
{
  typedef DualNumberBase<typename SymmetricOrType<T, T>::supertype, bool> supertype;
};

template <typename T, typename T2, typename D, bool reverseorder>
struct CompareTypes<DualNumberBase<T, D>,
                    T2,
                    reverseorder,
                    typename boostcopy::enable_if<BuiltinTraits<T2>>::type>
{
  typedef DualNumberBase<
      typename SymmetricCompareTypes<T, T2>::supertype,
      typename SymmetricCompareTypes<typename SymmetricCompareTypes<D, T2>::supertype,
                                     T>::supertype>
      supertype;
};

template <typename T, typename D, typename T2, typename D2>
struct CompareTypes<DualNumberBase<T, D>, DualNumberBase<T2, D2>>
{
  typedef DualNumberBase<
      typename SymmetricCompareTypes<T, T2>::supertype,
      typename SymmetricCompareTypes<typename SymmetricCompareTypes<T, T2>::supertype,
                                     typename SymmetricCompareTypes<D, D2>::supertype>::supertype>
      supertype;
};

template <typename T, typename D>
struct CompareTypes<DualNumberBase<T, D>, DualNumberBase<T, D>>
{
  typedef DualNumberBase<T, typename SymmetricCompareTypes<T, D>::supertype> supertype;
};

template <typename T, typename D>
inline D gradient(const DualNumberBase<T, D> & a);

} // namespace MetaPhysicL

namespace std
{

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DualNumberBase;

template <typename T, typename D>
inline bool isnan(const DualNumberBase<T, D> & a);

// Some forward declarations necessary for recursive DualNumberBases

#if __cplusplus >= 201103L

template <typename T, typename D>
inline DualNumberBase<T, D> cos(const DualNumberBase<T, D> & a);

template <typename T, typename D>
inline DualNumberBase<T, D> cos(DualNumberBase<T, D> && a);

template <typename T, typename D>
inline DualNumberBase<T, D> cosh(const DualNumberBase<T, D> & a);

template <typename T, typename D>
inline DualNumberBase<T, D> cosh(DualNumberBase<T, D> && a);

#else

template <typename T, typename D>
inline DualNumberBase<T, D> cos(DualNumberBase<T, D> a);

template <typename T, typename D>
inline DualNumberBase<T, D> cosh(DualNumberBase<T, D> a);

#endif

// Now just combined declaration/definitions

#if __cplusplus >= 201103L
#define DualNumber_decl_std_unary(funcname)                                                        \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(const DualNumberBase<T, D> & in);                           \
                                                                                                   \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(DualNumberBase<T, D> && in)

#else

#define DualNumber_decl_std_unary(funcname)                                                        \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(DualNumberBase<T, D> in)

#endif

DualNumber_decl_std_unary(sqrt);
DualNumber_decl_std_unary(exp);
DualNumber_decl_std_unary(log);
DualNumber_decl_std_unary(log10);
DualNumber_decl_std_unary(sin);
DualNumber_decl_std_unary(cos);
DualNumber_decl_std_unary(tan);
DualNumber_decl_std_unary(asin);
DualNumber_decl_std_unary(acos);
DualNumber_decl_std_unary(atan);
DualNumber_decl_std_unary(sinh);
DualNumber_decl_std_unary(cosh);
DualNumber_decl_std_unary(tanh);
DualNumber_decl_std_unary(abs);
DualNumber_decl_std_unary(fabs);
DualNumber_decl_std_unary(ceil);
DualNumber_decl_std_unary(floor);

#define DualNumber_decl_std_binary(funcname)                                                       \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline typename CompareTypes<DualNumberBase<T, D>, DualNumberBase<T2, D2>>::supertype funcname(  \
      const DualNumberBase<T, D> & a, const DualNumberBase<T2, D2> & b);                           \
                                                                                                   \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(const DualNumberBase<T, D> & a,                             \
                                       const DualNumberBase<T, D> & b);                            \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename CompareTypes<DualNumberBase<T2, D>, T, true>::supertype funcname(                \
      const T & a, const DualNumberBase<T2, D> & b);                                               \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename CompareTypes<DualNumberBase<T, D>, T2>::supertype funcname(                      \
      const DualNumberBase<T, D> & a, const T2 & b)

// if_else is necessary here to handle cases where a is negative but b
// is 0; we should have a contribution of 0 from those, not NaN.
DualNumber_decl_std_binary(pow);
DualNumber_decl_std_binary(atan2);
DualNumber_decl_std_binary(max);
DualNumber_decl_std_binary(min);
DualNumber_decl_std_binary(fmod);

template <typename T, typename D>
class numeric_limits<DualNumberBase<T, D>>
  : public MetaPhysicL::raw_numeric_limits<DualNumberBase<T, D>, T>
{
};

} // namespace std

#endif // METAPHYSICL_DUALNUMBER_DECL_H
