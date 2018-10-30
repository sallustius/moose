#ifndef MOOSEADWRAPPER_H
#define MOOSEADWRAPPER_H

#include "MooseTypes.h"
#include "MooseError.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"

#include <typeinfo>

template <typename T, bool declared_ad = false>
class MooseADWrapper
{
public:
  MooseADWrapper() : _val() {}

  typedef T DNType;

  const T & operator()(bool requested_by_user = true) const
  {
    if (requested_by_user)
      mooseWarning("Either type ",
                   typeid(T).name(),
                   " does not currently support automatic differentiation or the material property "
                   "was not declared as an AD property. Your Jacobian may be inaccurate");
    return _val;
  }

  T & operator()(bool requested_by_user = true)
  {
    if (requested_by_user)
      mooseWarning("Either type ",
                   typeid(T).name(),
                   " does not currently support automatic differentiation or the material property "
                   "was not declared as an AD property. Your Jacobian may be inaccurate");
    return _val;
  }

  void copyValueToDualNumber(const T & in) { _val = in; }
  void copyDualNumberToValue(T & in) { in = _val; }

private:
  T _val;
};

template <>
class MooseADWrapper<Real, true>
{
public:
  MooseADWrapper() : _val() {}

  typedef ADReal DNType;

  const ADReal & operator()(bool = true) const { return _val; }

  ADReal & operator()(bool = true) { return _val; }

  void copyValueToDualNumber(const Real & in) { _val.value() = in; }
  void copyDualNumberToValue(Real & in) { in = _val.value(); }

private:
  ADReal _val;
};

template <>
class MooseADWrapper<libMesh::VectorValue<Real>, true>
{
public:
  MooseADWrapper() : _val() {}

  typedef libMesh::VectorValue<ADReal> DNType;

  const libMesh::VectorValue<ADReal> & operator()(bool = true) const { return _val; }

  libMesh::VectorValue<ADReal> & operator()(bool = true) { return _val; }

  void copyValueToDualNumber(const VectorValue<Real> & in)
  {
    for (decltype(LIBMESH_DIM) i = 0; i < LIBMESH_DIM; ++i)
      _val(i).value() = in(i);
  }
  void copyDualNumberToValue(VectorValue<Real> & in)
  {
    for (decltype(LIBMESH_DIM) i = 0; i < LIBMESH_DIM; ++i)
      in(i) = _val(i).value();
  }

private:
  libMesh::VectorValue<ADReal> _val;
};

template <>
class MooseADWrapper<libMesh::TensorValue<Real>, true>
{
public:
  MooseADWrapper() : _val() {}

  typedef libMesh::TensorValue<ADReal> DNType;

  const libMesh::TensorValue<ADReal> & operator()(bool = true) const { return _val; }

  libMesh::TensorValue<ADReal> & operator()(bool = true) { return _val; }

  void copyValueToDualNumber(const TensorValue<Real> & in)
  {
    for (decltype(LIBMESH_DIM) i = 0; i < LIBMESH_DIM; ++i)
      for (decltype(LIBMESH_DIM) j = 0; j < LIBMESH_DIM; ++j)
        _val(i, j).value() = in(i, j);
  }
  void copyDualNumberToValue(TensorValue<Real> & in)
  {
    for (decltype(LIBMESH_DIM) i = 0; i < LIBMESH_DIM; ++i)
      for (decltype(LIBMESH_DIM) j = 0; j < LIBMESH_DIM; ++j)
        in(i, j) = _val(i, j).value();
  }

private:
  libMesh::TensorValue<ADReal> _val;
};

#endif // MOOSEADWRAPPER_H
