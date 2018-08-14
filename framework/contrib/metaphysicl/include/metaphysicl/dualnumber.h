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

#ifndef METAPHYSICL_DUALNUMBER_H
#define METAPHYSICL_DUALNUMBER_H

#include "metaphysicl/dualnumber_decl.h"

namespace MetaPhysicL
{

template <typename T, typename D>
inline T &
DualNumberBase<T, D>::value()
{
  return _val;
}

template <typename T, typename D>
inline const T &
DualNumberBase<T, D>::value() const
{
  return _val;
}

template <typename T, typename D>
inline D &
DualNumberBase<T, D>::derivatives()
{
  return _deriv;
}

template <typename T, typename D>
inline const D &
DualNumberBase<T, D>::derivatives() const
{
  return _deriv;
}

template <typename T, typename D>
inline bool
DualNumberBase<T, D>::boolean_test() const
{
  return _val;
}

template <typename T, typename D>
inline DualNumberBase<T, D>
DualNumberBase<T, D>::operator-() const
{
  return DualNumberBase<T, D>(-_val, -_deriv);
}

template <typename T, typename D>
inline DualNumberBase<T, D> DualNumberBase<T, D>::operator!() const
{
  return DualNumberBase<T, D>(!_val, !_deriv);
}

//
// Member function definitions
//

template <typename T, typename D>
inline DualNumberBase<T, D>::DualNumberBase() : _val(), _deriv()
{
}

template <typename T, typename D>
template <typename T2>
inline DualNumberBase<T, D>::DualNumberBase(const T2 & val)
  : _val(DualNumberConstructor<T, D>::value(val)), _deriv(DualNumberConstructor<T, D>::deriv(val))
{
}

template <typename T, typename D>
template <typename T2, typename D2>
inline DualNumberBase<T, D>::DualNumberBase(const T2 & val, const D2 & deriv)
  : _val(DualNumberConstructor<T, D>::value(val, deriv)),
    _deriv(DualNumberConstructor<T, D>::deriv(val, deriv))
{
}

/*
 * vector specialization definitions
 */
template <typename T, std::size_t N>
inline DualNumber<typename T::value_type, N>
DualNumber<T, N, typename std::enable_if<VectorTraits<T>::value>::type>::
operator()(const unsigned int i) const
{
  NumberArray<N, typename T::value_type> deriv;
  for (decltype(N) di = 0; i < N; ++di)
    deriv[di] = this->derivatives()[di](i);
  return {this->value()(i), deriv};
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<typename T::value_type, N> &
DualNumber<T, N, typename std::enable_if<VectorTraits<T>::value>::type>::
operator()(const unsigned int i)
{
  typename std::map<std::pair<unsigned int, unsigned int>,
                    DualNumberSurrogate<typename T::value_type, N>>::iterator it =
      _vector_dual_number_surrogates.find(i);
  if (it == _vector_dual_number_surrogates.end())
    _vector_dual_number_surrogates.emplace(i, {*this, i});
  return _vector_dual_number_surrogates[i];
}

template <typename T, std::size_t N>
inline void
DualNumber<T, N, typename std::enable_if<VectorTraits<T>::value>::type>::zero()
{
  zero_helper(*this);
}

/*
 * tensor specialization definitions
 */
template <typename T, std::size_t N>
inline auto
DualNumber<T, N, typename std::enable_if<TensorTraits<T>::value>::type>::row(
    const unsigned int r) const -> DualNumber<decltype(this->value().row(r)), N>
{
  NumberArray<N, decltype(this->value().row(r))> deriv;
  for (decltype(N) i = 0; i < N; ++i)
    deriv[i] = this->derivatives()[i].row(r);
  return {this->value().row(r), deriv};
}

template <typename T, std::size_t N>
inline DualNumber<typename T::value_type, N>
DualNumber<T, N, typename std::enable_if<TensorTraits<T>::value>::type>::tr() const
{
  NumberArray<N, typename T::value_type> deriv;
  for (decltype(N) i = 0; i < N; ++i)
    deriv[i] = this->derivatives()[i].tr();
  return {this->value().tr(), deriv};
}

template <typename T, std::size_t N>
inline DualNumber<typename T::value_type, N>
DualNumber<T, N, typename std::enable_if<TensorTraits<T>::value>::type>::
operator()(const unsigned int i, const unsigned int j) const
{
  NumberArray<N, typename T::value_type> deriv;
  for (decltype(N) di = 0; di < N; ++di)
    deriv[di] = this->derivatives()[di](i, j);
  return {this->value()(i, j), deriv};
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<typename T::value_type, N> &
DualNumber<T, N, typename std::enable_if<TensorTraits<T>::value>::type>::
operator()(const unsigned int i, const unsigned int j)
{
  std::pair<unsigned, unsigned> indices(i, j);
  typename std::map<std::pair<unsigned int, unsigned int>,
                    DualNumberSurrogate<typename T::value_type, N>>::iterator it =
      _tensor_dual_number_surrogates.find(indices);
  if (it == _tensor_dual_number_surrogates.end())
    _tensor_dual_number_surrogates.emplace(indices, {*this, i, j});
  return _tensor_dual_number_surrogates[indices];
}

template <typename T, std::size_t N>
inline void
DualNumber<T, N, typename std::enable_if<TensorTraits<T>::value>::type>::zero()
{
  zero_helper(*this);
}

/*
 * RankFourTensor specialization definitions
 */
template <typename T, std::size_t N>
inline DualNumber<Real, N>
DualNumber<T, N, typename std::enable_if<RankFourTraits<T>::value>::type>::operator()(
    const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l) const
{
  NumberArray<N, Real> deriv;
  for (decltype(N) di = 0; di < N; ++di)
    deriv[di] = this->derivatives()[di](i, j, k, l);
  return {this->value()(i, j, k, l), deriv};
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<Real, N> &
DualNumber<T, N, typename std::enable_if<RankFourTraits<T>::value>::type>::
operator()(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l)
{
  std::tuple<unsigned, unsigned, unsigned, unsigned> indices(i, j, k, l);
  typename std::map<std::tuple<unsigned, unsigned, unsigned, unsigned>,
                    DualNumberSurrogate<Real, N>>::iterator it =
      _rank4_dual_number_surrogates.find(indices);
  if (it == _rank4_dual_number_surrogates.end())
    _rank4_dual_number_surrogates.emplace(indices, {*this, i, j, k, l});
  return _rank4_dual_number_surrogates[indices];
}

// helper function
template <typename T>
void
zero_helper(T & dn)
{
  dn.value().zero();
  auto size = dn.derivatives().size();
  for (decltype(size) i = 0; i < size; ++i)
    dn.derivatives()[i].zero();
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<T, N>::DualNumberSurrogate(DualNumber<T, N> & dn) : value(dn.value())
{
  for (decltype(N) i = 0; i < N; ++i)
    derivatives[i] = &dn.derivatives()[i];
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<T, N>::DualNumberSurrogate(DualNumber<T, N> && dn) : value(dn.value())
{
  for (decltype(N) i = 0; i < N; ++i)
    derivatives[i] = &dn.derivatives()[i];
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<T, N>::DualNumberSurrogate(const T & n) : value(n)
{
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<T, N>::DualNumberSurrogate(const DualNumberSurrogate<T, N> & dns)
  : value(dns.value), derivatives(dns.derivatives)
{
}

template <typename T, std::size_t N>
inline DualNumberSurrogate<T, N> &
DualNumberSurrogate<T, N>::operator=(const DualNumberSurrogate<T, N> & dns)
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

  // FIXME: these operators currently do automatic type promotion when
  // encountering DualNumbers of differing levels of recursion and
  // differentiability.  But what we really want is automatic type
  // *demotion*, to avoid pretending we have accurate derivatives which
  // we don't have.  If we could do that right then some potential
  // subtle run-time user errors would turn into compile-time user
  // errors.

#define DualNumber_op(opname, functorname, dn_first_calc, dn_second_calc, dualcalc)                \
  template <typename T, typename D>                                                                \
  template <typename T2>                                                                           \
  inline auto DualNumberBase<T, D>::operator opname##=(const T2 & in)->decltype(*this) &           \
  {                                                                                                \
    const auto & a = *this;                                                                        \
    const auto & b = in;                                                                           \
    this->derivatives() = dn_first_calc;                                                           \
    this->value() opname## = b;                                                                    \
    return *this;                                                                                  \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D>                                                                \
  template <typename T2, typename D2>                                                              \
  inline auto DualNumberBase<T, D>::operator opname##=(const DualNumberBase<T2, D2> & in)          \
      ->decltype(*this) &                                                                          \
  {                                                                                                \
    const auto & a = *this;                                                                        \
    const auto & b = in;                                                                           \
    this->derivatives() = dualcalc;                                                                \
    this->value() opname## = b.value();                                                            \
    return *this;                                                                                  \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D>                                                                \
  template <typename T2, size_t N>                                                                 \
  inline auto DualNumberBase<T, D>::operator opname##=(const DualNumber<T2, N> & in)               \
      ->decltype(*this) &                                                                          \
  {                                                                                                \
    const auto & a = *this;                                                                        \
    const auto & b = in;                                                                           \
    this->derivatives() = dualcalc;                                                                \
    this->value() opname## = b.value();                                                            \
    return *this;                                                                                  \
  }                                                                                                \
                                                                                                   \
  template <typename T, size_t N, typename T2>                                                     \
  inline auto operator opname(const DualNumber<T, N> & a, const DualNumber<T2, N> & b)             \
      ->DualNumber<decltype(a.value() opname b.value()), N>                                        \
  {                                                                                                \
    auto value = a.value() opname b.value();                                                       \
    auto derivatives = dualcalc;                                                                   \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, size_t N>                                                     \
  inline auto operator opname(const T & a, const DualNumber<T2, N> & b)                            \
      ->DualNumber<decltype(a opname b.value()), N>                                                \
  {                                                                                                \
    auto value = a opname b.value();                                                               \
    auto derivatives = dn_second_calc;                                                             \
    return {value, derivatives};                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, size_t N, typename T2>                                                     \
  inline auto operator opname(const DualNumber<T, N> & a, const T2 & b)                            \
      ->DualNumber<decltype(a.value() opname b), N>                                                \
  {                                                                                                \
    auto value = a.value() opname b;                                                               \
    auto derivatives = dn_first_calc;                                                              \
    return {value, derivatives};                                                                   \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

DualNumber_op(+, Plus, a.derivatives(), b.derivatives(), a.derivatives() + b.derivatives());

DualNumber_op(-, Minus, a.derivatives(), -b.derivatives, a.derivatives() - b.derivatives());

DualNumber_op(*,
              Multiplies,
              a.derivatives() * b,
              a * b.derivatives(),
              a.value() * b.derivatives() + a.derivatives() * b.value());

DualNumber_op(/,
              Divides,
              a.derivatives() / b,
              -a * b.derivatives() / (b.value() * b.value()),
              (b.value() * a.derivatives() - b.derivatives() * a.value()) /
                  (b.value() * b.value()));

#define DualNumber_compare(opname)                                                                 \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline bool operator opname(const DualNumberBase<T, D> & a, const DualNumberBase<T2, D2> & b)    \
  {                                                                                                \
    return (a.value() opname b.value());                                                           \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, typename D2>                                                  \
  inline typename boostcopy::                                                                      \
      enable_if_class<typename CompareTypes<DualNumberBase<T2, D2>, T>::supertype, bool>::type     \
      operator opname(const T & a, const DualNumberBase<T2, D2> & b)                               \
  {                                                                                                \
    return (a opname b.value());                                                                   \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename boostcopy::                                                                      \
      enable_if_class<typename CompareTypes<DualNumberBase<T, D>, T2>::supertype, bool>::type      \
      operator opname(const DualNumberBase<T, D> & a, const T2 & b)                                \
  {                                                                                                \
    return (a.value() opname b);                                                                   \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

DualNumber_compare(>);
DualNumber_compare(>=);
DualNumber_compare(<);
DualNumber_compare(<=);
DualNumber_compare(==);
DualNumber_compare(!=);
DualNumber_compare(&&);
DualNumber_compare(||);

template <typename T, typename D>
inline std::ostream &
operator<<(std::ostream & output, const DualNumberBase<T, D> & a)
{
  return output << '(' << a.value() << ',' << a.derivatives() << ')';
}

template <typename T, typename D>
inline D
gradient(const DualNumberBase<T, D> & a)
{
  return a.derivatives();
}

} // namespace MetaPhysicL

namespace std
{

using MetaPhysicL::CompareTypes;
using MetaPhysicL::DualNumberBase;

template <typename T, typename D>
inline bool
isnan(const DualNumberBase<T, D> & a)
{
  using std::isnan;
  return isnan(a.value());
}

#if __cplusplus >= 201103L
#define DualNumber_std_unary(funcname, derivative, precalc)                                        \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(const DualNumberBase<T, D> & in)                            \
  {                                                                                                \
    DualNumberBase<T, D> returnval = in;                                                           \
    T funcval = std::funcname(in.value());                                                         \
    precalc;                                                                                       \
    returnval.derivatives() *= derivative;                                                         \
    returnval.value() = funcval;                                                                   \
    return returnval;                                                                              \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(DualNumberBase<T, D> && in)                                 \
  {                                                                                                \
    T funcval = std::funcname(in.value());                                                         \
    precalc;                                                                                       \
    in.derivatives() *= derivative;                                                                \
    in.value() = funcval;                                                                          \
    return in;                                                                                     \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

#else

#define DualNumber_std_unary(funcname, derivative, precalc)                                        \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(DualNumberBase<T, D> in)                                    \
  {                                                                                                \
    T funcval = std::funcname(in.value());                                                         \
    precalc;                                                                                       \
    in.derivatives() *= derivative;                                                                \
    in.value() = funcval;                                                                          \
    return in;                                                                                     \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

#endif

DualNumber_std_unary(sqrt, 1 / (2 * funcval), );
DualNumber_std_unary(exp, funcval, );
DualNumber_std_unary(log, 1 / in.value(), );
DualNumber_std_unary(log10, 1 / in.value() * (1 / std::log(T(10.))), );
DualNumber_std_unary(sin, std::cos(in.value()), );
DualNumber_std_unary(cos, -std::sin(in.value()), );
DualNumber_std_unary(tan, sec_in * sec_in, T sec_in = 1 / std::cos(in.value()));
DualNumber_std_unary(asin, 1 / std::sqrt(1 - in.value() * in.value()), );
DualNumber_std_unary(acos, -1 / std::sqrt(1 - in.value() * in.value()), );
DualNumber_std_unary(atan, 1 / (1 + in.value() * in.value()), );
DualNumber_std_unary(sinh, std::cosh(in.value()), );
DualNumber_std_unary(cosh, std::sinh(in.value()), );
DualNumber_std_unary(tanh, sech_in * sech_in, T sech_in = 1 / std::cosh(in.value()));
DualNumber_std_unary(abs, (in.value() > 0) - (in.value() < 0), );  // std < and > return 0 or 1
DualNumber_std_unary(fabs, (in.value() > 0) - (in.value() < 0), ); // std < and > return 0 or 1
DualNumber_std_unary(ceil, 0, );
DualNumber_std_unary(floor, 0, );

#define DualNumber_std_binary(funcname, derivative)                                                \
  template <typename T, typename D, typename T2, typename D2>                                      \
  inline typename CompareTypes<DualNumberBase<T, D>, DualNumberBase<T2, D2>>::supertype funcname(  \
      const DualNumberBase<T, D> & a, const DualNumberBase<T2, D2> & b)                            \
  {                                                                                                \
    typedef typename CompareTypes<T, T2>::supertype TS;                                            \
    typedef typename CompareTypes<DualNumberBase<T, D>, DualNumberBase<T2, D2>>::supertype type;   \
                                                                                                   \
    TS funcval = std::funcname(a.value(), b.value());                                              \
    return type(funcval, derivative);                                                              \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename D>                                                                \
  inline DualNumberBase<T, D> funcname(const DualNumberBase<T, D> & a,                             \
                                       const DualNumberBase<T, D> & b)                             \
  {                                                                                                \
    T funcval = std::funcname(a.value(), b.value());                                               \
    return DualNumberBase<T, D>(funcval, derivative);                                              \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename CompareTypes<DualNumberBase<T2, D>, T, true>::supertype funcname(                \
      const T & a, const DualNumberBase<T2, D> & b)                                                \
  {                                                                                                \
    typedef typename CompareTypes<DualNumberBase<T2, D>, T, true>::supertype type;                 \
    type newa(a);                                                                                  \
    return std::funcname(newa, b);                                                                 \
  }                                                                                                \
                                                                                                   \
  template <typename T, typename T2, typename D>                                                   \
  inline typename CompareTypes<DualNumberBase<T, D>, T2>::supertype funcname(                      \
      const DualNumberBase<T, D> & a, const T2 & b)                                                \
  {                                                                                                \
    typedef typename CompareTypes<DualNumberBase<T, D>, T2>::supertype type;                       \
    type newb(b);                                                                                  \
    return std::funcname(a, newb);                                                                 \
  }                                                                                                \
  void ANONYMOUS_FUNCTION()

// if_else is necessary here to handle cases where a is negative but b
// is 0; we should have a contribution of 0 from those, not NaN.
DualNumber_std_binary(pow,
                      funcval *(b.value() * a.derivatives() / a.value() +
                                MetaPhysicL::if_else(b.derivatives(),
                                                     b.derivatives() * std::log(a.value()),
                                                     b.derivatives())));
DualNumber_std_binary(atan2,
                      (b.value() * a.derivatives() - a.value() * b.derivatives()) /
                          (b.value() * b.value() + a.value() * a.value()));
DualNumber_std_binary(max, (a.value() > b.value()) ? a.derivatives() : b.derivatives());
DualNumber_std_binary(min, (a.value() > b.value()) ? b.derivatives() : a.derivatives());
DualNumber_std_binary(fmod, a.derivatives());

} // namespace std

#endif // METAPHYSICL_DUALNUMBER_H
