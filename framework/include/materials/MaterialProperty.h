//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef MATERIALPROPERTY_H
#define MATERIALPROPERTY_H

#include <vector>

#include "MooseArray.h"
#include "DataIO.h"
#include "MooseADWrapper.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

class PropertyValue;

/**
 * Abstract definition of a property value.
 */
class PropertyValue
{
public:
  virtual ~PropertyValue(){};

  /**
   * String identifying the type of parameter stored.
   */
  virtual std::string type() = 0;

  /**
   * Clone this value.  Useful in copy-construction.
   */
  virtual PropertyValue * init(int size) = 0;

  virtual unsigned int size() const = 0;

  /**
   * Resizes the property to the size n
   */
  virtual void resize(int n) = 0;

  virtual void swap(PropertyValue * rhs) = 0;

  /**
   * Copy the value of a Property from one specific to a specific qp in this Property.
   *
   * @param to_qp The quadrature point in _this_ Property that you want to copy to.
   * @param rhs The Property you want to copy _from_.
   * @param from_qp The quadrature point in rhs you want to copy _from_.
   */
  virtual void
  qpCopy(const unsigned int to_qp, PropertyValue * rhs, const unsigned int from_qp) = 0;

  // save/restore in a file
  virtual void store(std::ostream & stream) = 0;
  virtual void load(std::istream & stream) = 0;

  /**
   * copy the Real version to the DualNumber<Real> version of the material property for the
   * specified quadrature point
   */
  virtual void copyValueToDualNumber(const unsigned int i) = 0;

  /**
   * copy the value portion (not the derivatives) of the DualNumber<Real> version of the material
   * property to the Real version for the specified quadrature point
   */
  virtual void copyDualNumberToValue(const unsigned int i) = 0;
};

template <>
inline void
dataStore(std::ostream & stream, PropertyValue *& p, void * /*context*/)
{
  p->store(stream);
}

template <>
inline void
dataLoad(std::istream & stream, PropertyValue *& p, void * /*context*/)
{
  p->load(stream);
}

/**
 * Concrete definition of a parameter value
 * for a specified type.
 */
template <typename T>
class MaterialProperty : public PropertyValue
{
public:
  /// Explicitly declare a public constructor because we made the copy constructor private
  MaterialProperty() : PropertyValue()
  { /* */
  }

  virtual ~MaterialProperty() { _value.release(); }

  /**
   * @returns a read-only reference to the parameter value.
   */
  const MooseArray<T> & get() const { return _value; }

  /**
   * @returns a writable reference to the parameter value.
   */
  MooseArray<T> & set() { return _value; }

  /**
   * String identifying the type of parameter stored.
   */
  virtual std::string type() override;

  /**
   * Clone this value.  Useful in copy-construction.
   */
  virtual PropertyValue * init(int size) override;

  /**
   * Resizes the property to the size n
   */
  virtual void resize(int n) override;

  virtual unsigned int size() const override { return _value.size(); }

  /**
   * Get element i out of the array as a writeable reference.
   */
  T & operator[](const unsigned int i) { return _value[i]; }

  /**
   * Get element i out of the array as a ready-only reference.
   */
  const T & operator[](const unsigned int i) const { return _value[i]; }

  /**
   *
   */
  virtual void swap(PropertyValue * rhs) override;

  /**
   * Copy the value of a Property from one specific to a specific qp in this Property.
   *
   * @param to_qp The quadrature point in _this_ Property that you want to copy to.
   * @param rhs The Property you want to copy _from_.
   * @param from_qp The quadrature point in rhs you want to copy _from_.
   */
  virtual void
  qpCopy(const unsigned int to_qp, PropertyValue * rhs, const unsigned int from_qp) override;

  /**
   * Store the property into a binary stream
   */
  virtual void store(std::ostream & stream) override;

  /**
   * Load the property from a binary stream
   */
  virtual void load(std::istream & stream) override;

  void copyValueToDualNumber(const unsigned int /*i*/) override {}

  void copyDualNumberToValue(const unsigned int /*i*/) override {}

private:
  /// private copy constructor to avoid shallow copying of material properties
  MaterialProperty(const MaterialProperty<T> & /*src*/)
  {
    mooseError("Material properties must be assigned to references (missing '&')");
  }

  /// private assignment operator to avoid shallow copying of material properties
  MaterialProperty<T> & operator=(const MaterialProperty<T> & /*rhs*/)
  {
    mooseError("Material properties must be assigned to references (missing '&')");
  }

protected:
  /// Stored parameter value.
  MooseArray<T> _value;
};

// ------------------------------------------------------------
// Material::Property<> class inline methods
template <typename T>
inline std::string
MaterialProperty<T>::type()
{
  return typeid(T).name();
}

template <typename T>
inline PropertyValue *
MaterialProperty<T>::init(int size)
{
  MaterialProperty<T> * copy = new MaterialProperty<T>;
  copy->_value.resize(size, T{});
  return copy;
}

template <typename T>
inline void
MaterialProperty<T>::resize(int n)
{
  _value.resize(n);
}

template <typename T>
inline void
MaterialProperty<T>::swap(PropertyValue * rhs)
{
  mooseAssert(rhs != NULL, "Assigning NULL?");
  _value.swap(cast_ptr<MaterialProperty<T> *>(rhs)->_value);
}

template <typename T>
inline void
MaterialProperty<T>::qpCopy(const unsigned int to_qp,
                            PropertyValue * rhs,
                            const unsigned int from_qp)
{
  mooseAssert(rhs != NULL, "Assigning NULL?");
  _value[to_qp] = cast_ptr<const MaterialProperty<T> *>(rhs)->_value[from_qp];
}

template <typename T>
inline void
MaterialProperty<T>::store(std::ostream & stream)
{
  for (unsigned int i = 0; i < size(); i++)
    storeHelper(stream, _value[i], NULL);
}

template <typename T>
inline void
MaterialProperty<T>::load(std::istream & stream)
{
  for (unsigned int i = 0; i < size(); i++)
    loadHelper(stream, _value[i], NULL);
}

template <typename T>
class ADMaterialPropertyObject : public MaterialProperty<T>
{
public:
  ADMaterialPropertyObject() : MaterialProperty<T>() {}

  virtual ~ADMaterialPropertyObject() { _dual_number.release(); }

  /**
   * @returns a read-only reference to the parameter dual number.
   */
  const MooseArray<MooseADWrapper<T>> & get() const { return _dual_number; }

  /**
   * @returns a writable reference to the parameter dual number.
   */
  MooseArray<MooseADWrapper<T>> & set() { return _dual_number; }

  virtual PropertyValue * init(int size) override;
  virtual void resize(int n) override;
  virtual void swap(PropertyValue * rhs) override;

  virtual void
  qpCopy(const unsigned int to_qp, PropertyValue * rhs, const unsigned int from_qp) override;
  virtual void store(std::ostream & stream) override;
  virtual void load(std::istream & stream) override;

  /**
   * Get element i out of the array as a writeable reference.
   */
  typename MooseADWrapper<T>::DNType & operator[](const unsigned int i)
  {
    return _dual_number[i]();
  }
  /**
   * Get element i out of the array as a read-only reference.
   */
  const typename MooseADWrapper<T>::DNType & operator[](const unsigned int i) const
  {
    return _dual_number[i]();
  }

  /**
   * copy the Real version to the DualNumber<Real> version of the material property for the
   * specified quadrature point
   */
  void copyValueToDualNumber(const unsigned int i) override
  {
    _dual_number[i].copyValueToDualNumber(this->_value[i]);
  }

  /**
   * copy the value portion (not the derivatives) of the DualNumber<Real> version of the material
   * property to the Real version for the specified quadrature point
   */
  void copyDualNumberToValue(const unsigned int i) override
  {
    _dual_number[i].copyDualNumberToValue(this->_value[i]);
  }

protected:
  /// Stored dual number
  MooseArray<MooseADWrapper<T>> _dual_number;
};

template <typename T>
inline PropertyValue *
ADMaterialPropertyObject<T>::init(int size)
{
  ADMaterialPropertyObject<T> * copy = new ADMaterialPropertyObject<T>;
  copy->_value.resize(size, T{});
  copy->_dual_number.resize(size, MooseADWrapper<T>{});
  return copy;
}

template <typename T>
inline void
ADMaterialPropertyObject<T>::resize(int n)
{
  _dual_number.resize(n);
  this->_value.resize(n);
}

template <typename T>
inline void
ADMaterialPropertyObject<T>::swap(PropertyValue * rhs)
{
  mooseAssert(rhs != NULL, "Assigning NULL?");
  this->_value.swap(cast_ptr<ADMaterialPropertyObject<T> *>(rhs)->_value);
  _dual_number.swap(cast_ptr<ADMaterialPropertyObject<T> *>(rhs)->_dual_number);
}

template <typename T>
inline void
ADMaterialPropertyObject<T>::qpCopy(const unsigned int to_qp,
                                    PropertyValue * rhs,
                                    const unsigned int from_qp)
{
  mooseAssert(rhs != NULL, "Assigning NULL?");
  this->_value[to_qp] = cast_ptr<const ADMaterialPropertyObject<T> *>(rhs)->_value[from_qp];
  _dual_number[to_qp] = cast_ptr<const ADMaterialPropertyObject<T> *>(rhs)->_dual_number[from_qp];
}

template <typename T>
inline void
ADMaterialPropertyObject<T>::store(std::ostream & stream)
{
  for (unsigned int i = 0; i < this->size(); i++)
  {
    storeHelper(stream, this->_value[i], nullptr);
    storeHelper(stream, _dual_number[i], nullptr);
  }
}

template <typename T>
inline void
ADMaterialPropertyObject<T>::load(std::istream & stream)
{
  for (unsigned int i = 0; i < this->size(); i++)
  {
    loadHelper(stream, this->_value[i], nullptr);
    loadHelper(stream, _dual_number[i], nullptr);
  }
}

/**
 * Container for storing material properties
 */
class MaterialProperties : public std::vector<PropertyValue *>
{
public:
  MaterialProperties() {}

  virtual ~MaterialProperties() {}

  /**
   * Parameter map iterator.
   */
  typedef std::vector<PropertyValue *>::iterator iterator;

  /**
   * Constant parameter map iterator.
   */
  typedef std::vector<PropertyValue *>::const_iterator const_iterator;

  /**
   * Deallocates the memory
   */
  void destroy()
  {
    for (iterator k = begin(); k != end(); ++k)
      delete *k;
  }

  /**
   * Resize items in this array, i.e. the number of values needed in PropertyValue array
   * @param n_qpoints The number of values needed to store (equals the the number of quadrature
   * points per mesh element)
   */
  void resizeItems(unsigned int n_qpoints)
  {
    for (iterator k = begin(); k != end(); ++k)
      if (*k != NULL)
        (*k)->resize(n_qpoints);
  }
};

template <>
inline void
dataStore(std::ostream & stream, MaterialProperties & v, void * context)
{
  // Cast this to a vector so we can just piggy back on the vector store capability
  std::vector<PropertyValue *> & mat_props = static_cast<std::vector<PropertyValue *> &>(v);

  storeHelper(stream, mat_props, context);
}

template <>
inline void
dataLoad(std::istream & stream, MaterialProperties & v, void * context)
{
  // Cast this to a vector so we can just piggy back on the vector store capability
  std::vector<PropertyValue *> & mat_props = static_cast<std::vector<PropertyValue *> &>(v);

  loadHelper(stream, mat_props, context);
}

#endif
