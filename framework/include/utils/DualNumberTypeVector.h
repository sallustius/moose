#ifndef DUALNUMBERTYPEVECTOR_H
#define DUALNUMBERTYPEVECTOR_H

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_METAPHYSICL

#include "libmesh/type_vector.h"
#include "metaphysicl/compare_types.h"

template <typename T, bool reverseorder>
struct MetaPhysicL::PlusType<libMesh::TypeVector<T>, libMesh::TypeVector<T>, reverseorder>
{
  typedef libMesh::TypeVector<T> supertype;
};

template <typename T, bool reverseorder>
struct MetaPhysicL::MinusType<libMesh::TypeVector<T>, libMesh::TypeVector<T>, reverseorder>
{
  typedef libMesh::TypeVector<T> supertype;
};

template <typename T, bool reverseorder>
struct MetaPhysicL::MultipliesType<libMesh::TypeVector<T>, libMesh::TypeVector<T>, reverseorder>
{
  typedef Real supertype;
};

template <typename T, bool reverseorder>
struct MetaPhysicL::MultipliesType<Real, libMesh::TypeVector<T>, reverseorder>
{
  typedef libMesh::TypeVector<T> supertype;
};

template <typename T>
struct MetaPhysicL::DividesType<libMesh::TypeVector<T>, Real>
{
  typedef libMesh::TypeVector<T> supertype;
};

#endif
#endif
