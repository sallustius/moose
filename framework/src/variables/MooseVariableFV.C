//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseVariableFV.h"
#include <typeinfo>
#include "TimeIntegrator.h"
#include "NonlinearSystemBase.h"
#include "DisplacedSystem.h"
#include "Assembly.h"

#include "libmesh/numeric_vector.h"

registerMooseObject("MooseApp", MooseVariableFVReal);

template <typename OutputType>
InputParameters
MooseVariableFV<OutputType>::validParams()
{
  return MooseVariableField<OutputType>::validParams();
}

template <typename OutputType>
MooseVariableFV<OutputType>::MooseVariableFV(const InputParameters & parameters)
  : MooseVariableField<OutputType>(parameters),
    _solution(_sys.currentSolution()),
    _phi(&_assembly.template fePhi<OutputShape>(FEType(CONSTANT, MONOMIAL)))
{
  _element_data = libmesh_make_unique<MooseVariableDataFV<OutputType>>(
      *this, _sys, _tid, Moose::ElementType::Element, this->_assembly.elem());
  _neighbor_data = libmesh_make_unique<MooseVariableDataFV<OutputType>>(
      *this, _sys, _tid, Moose::ElementType::Neighbor, this->_assembly.neighbor());
}

template <typename OutputType>
const std::set<SubdomainID> &
MooseVariableFV<OutputType>::activeSubdomains() const
{
  return this->_sys.system().variable(_var_num).active_subdomains();
}

template <typename OutputType>
Moose::VarFieldType
MooseVariableFV<OutputType>::fieldType() const
{
  if (std::is_same<OutputType, Real>::value)
    return Moose::VarFieldType::VAR_FIELD_STANDARD;
  else if (std::is_same<OutputType, RealVectorValue>::value)
    return Moose::VarFieldType::VAR_FIELD_VECTOR;
  else if (std::is_same<OutputType, RealEigenVector>::value)
    return Moose::VarFieldType::VAR_FIELD_ARRAY;
  else
    mooseError("Unknown variable field type");
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::activeOnSubdomain(SubdomainID subdomain) const
{
  return this->_sys.system().variable(_var_num).active_on_subdomain(subdomain);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::clearDofIndices()
{
  _element_data->clearDofIndices();
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValue(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Current, idx);
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValueOld(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Old, idx);
}

template <typename OutputType>
typename MooseVariableFV<OutputType>::OutputData
MooseVariableFV<OutputType>::getElementalValueOlder(const Elem * elem, unsigned int idx) const
{
  return _element_data->getElementalValue(elem, Moose::Older, idx);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::insert(NumericVector<Number> & residual)
{
  _element_data->insert(residual);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::add(NumericVector<Number> & /*residual*/)
{
  mooseError("add not supported for FV variables");
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValues() const
{
  return _element_data->dofValues();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOld() const
{
  return _element_data->dofValuesOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOlder() const
{
  return _element_data->dofValuesOlder();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesPreviousNL() const
{
  return _element_data->dofValuesPreviousNL();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesNeighbor() const
{
  return _neighbor_data->dofValues();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOldNeighbor() const
{
  return _neighbor_data->dofValuesOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesOlderNeighbor() const
{
  return _neighbor_data->dofValuesOlder();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesPreviousNLNeighbor() const
{
  return _neighbor_data->dofValuesPreviousNL();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDot() const
{
  return _element_data->dofValuesDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDot() const
{
  return _element_data->dofValuesDotDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotOld() const
{
  return _element_data->dofValuesDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotOld() const
{
  return _element_data->dofValuesDotDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotNeighbor() const
{
  return _neighbor_data->dofValuesDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotNeighbor() const
{
  return _neighbor_data->dofValuesDotDot();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotOldNeighbor() const
{
  return _neighbor_data->dofValuesDotOld();
}

template <typename OutputType>
const typename MooseVariableFV<OutputType>::DoFValue &
MooseVariableFV<OutputType>::dofValuesDotDotOldNeighbor() const
{
  return _neighbor_data->dofValuesDotDotOld();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDu() const
{
  return _element_data->dofValuesDuDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDotDu() const
{
  return _element_data->dofValuesDuDotDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDuNeighbor() const
{
  return _neighbor_data->dofValuesDuDotDu();
}

template <typename OutputType>
const MooseArray<Number> &
MooseVariableFV<OutputType>::dofValuesDuDotDotDuNeighbor() const
{
  return _neighbor_data->dofValuesDuDotDotDu();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::prepareIC()
{
  _element_data->prepareIC();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeElemValues()
{
  _element_data->setGeometry(Moose::Volume);
  _element_data->computeValues();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::computeFaceValues(const FaceInfo & fi)
{
  _element_data->setGeometry(Moose::Face);
  _neighbor_data->setGeometry(Moose::Face);

  auto facetype = fi.faceType(_var_name);
  if (facetype == FaceInfo::VarFaceNeighbors::NEITHER)
    return;
  else if (facetype == FaceInfo::VarFaceNeighbors::BOTH)
  {
    _element_data->computeValuesFace(fi);
    _neighbor_data->computeValuesFace(fi);
  }
  else if (facetype == FaceInfo::VarFaceNeighbors::ELEM)
  {
    _element_data->computeValuesFace(fi);
    _neighbor_data->computeGhostValuesFace(fi, *_element_data);
  }
  else if (facetype == FaceInfo::VarFaceNeighbors::NEIGHBOR)
  {
    _neighbor_data->computeValuesFace(fi);
    _element_data->computeGhostValuesFace(fi, *_neighbor_data);
  }
  else
    mooseError("robert wrote broken MooseVariableFV code");
}

template <typename OutputType>
OutputType
MooseVariableFV<OutputType>::getValue(const Elem * elem) const
{
  std::vector<dof_id_type> dof_indices;
  this->_dof_map.dof_indices(elem, dof_indices, _var_num);
  mooseAssert(dof_indices.size() == 1, "Wrong size for dof indices");
  OutputType value = (*this->_sys.currentSolution())(dof_indices[0]);
  return value;
}

template <typename OutputType>
typename OutputTools<OutputType>::OutputGradient
MooseVariableFV<OutputType>::getGradient(const Elem * /*elem*/) const
{
  return {};
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setNodalValue(const OutputType & /*value*/, unsigned int /*idx*/)
{
  mooseError("FV variables do not support setNodalValue");
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setDofValue(const OutputData & value, unsigned int index)
{
  _element_data->setDofValue(value, index);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setElementalValue(const OutputType & value)
{
  _element_data->setElementalValue(value);
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::setDofValues(const DenseVector<OutputData> & values)
{
  _element_data->setDofValues(values);
}

template <typename OutputType>
bool
MooseVariableFV<OutputType>::isVector() const
{
  return std::is_same<OutputType, RealVectorValue>::value;
}

#ifdef MOOSE_GLOBAL_AD_INDEXING

template <typename OutputType>
ADReal
MooseVariableFV<OutputType>::getElemValue(const Elem * const elem) const
{
  std::vector<dof_id_type> dof_indices;
  _dof_map.dof_indices(elem, dof_indices, _var_num);

  mooseAssert(
      dof_indices.size() == 1,
      "There should only be one dof-index for a constant monomial variable on any given element");

  dof_id_type index = dof_indices[0];

  ADReal value = (*_solution)(index);

  if (ADReal::do_derivatives)
    Moose::derivInsert(value.derivatives(), index, 1.);

  return value;
}

template <typename OutputType>
const VectorValue<ADReal> &
MooseVariableFV<OutputType>::adGradSln(const Elem * const elem) const
{
  auto it = _elem_to_grad.find(elem);

  if (it != _elem_to_grad.end())
    return it->second;

  // Returns a pair with the first being an iterator pointing to the key-value pair and the second a
  // boolean denoting whether a new insertion took place
  auto emplace_ret = _elem_to_grad.emplace(elem, 0);

  mooseAssert(emplace_ret.second, "We should have inserted a new key-value pair");

  VectorValue<ADReal> & grad = emplace_ret.first->second;

  bool volume_set = false;
  Real volume = 0;

  ADReal elem_value = getElemValue(elem);

  for (const auto side : elem->side_index_range())
  {
    const Elem * const neighbor = elem->neighbor_ptr(side);

    bool elem_has_info = neighbor ? (elem->id() < neighbor->id()) : true;

    const FaceInfo * const fi = elem_has_info
                                    ? _mesh.faceInfo(elem, side)
                                    : _mesh.faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));

    mooseAssert(fi, "We should have found a FaceInfo");

    // Below uses notation from Moukalled's Finite Volume Method in Computational Fluid Dynamics,
    // Green-Gauss computation of the gradient

    auto face_value_functor = [&]() {
      if (neighbor)
      {
        ADReal neighbor_value = getElemValue(neighbor);

        const Real g_C = fi->gC();

        return g_C * elem_value + (1. - g_C) * neighbor_value;
      }
      else
      {
        // If we don't have a neighbor, then we're along a boundary, and we may have a DirichletBC

        std::vector<FVDirichletBC *> bcs;

        // TODO: this query probably (maybe?)needs to also filter based on the
        // active tags - these currently live in the flux thread loop object and I'm
        // not sure how best to get them here.
        _subproblem.getMooseApp()
            .theWarehouse()
            .query()
            .template condition<AttribSystem>("FVDirichletBC")
            .template condition<AttribThread>(_tid)
            .template condition<AttribBoundaries>(fi->boundaryIDs())
            .template condition<AttribVar>(_var_num)
            .queryInto(bcs);
        mooseAssert(bcs.size() <= 1, "cannot have multiple dirichlet BCs on the same boundary");

        bool has_dirichlet_bc = bcs.size() > 0;

        if (has_dirichlet_bc)
        {
          const FVDirichletBC & bc = *bcs[0];

          return ADReal(bc.boundaryValue(*fi));
        }
        else
          // No DirichletBC so we'll implicitly apply a zero gradient condition and assume that the
          // face value is equivalent to the element value
          return elem_value;
      }
    };

    const Point elem_normal = elem_has_info ? fi->normal() : Point(-fi->normal());
    const Point surface_vector = elem_normal * fi->faceArea();

    grad += face_value_functor() * surface_vector;

    if (!volume_set)
    {
      volume = elem_has_info ? fi->elemVolume() : fi->neighborVolume();
      volume_set = true;
    }
  }

  mooseAssert(volume_set && volume > 0, "We should have set the volume");

  grad /= volume;

  return grad;
}

template <typename OutputType>
const VectorValue<ADReal> &
MooseVariableFV<OutputType>::uncorrectedAdGradSln(const FaceInfo & fi) const
{
  auto it = _face_to_unc_grad.find(&fi);

  if (it != _face_to_unc_grad.end())
    return it->second;

  const VectorValue<ADReal> & elem_grad = adGradSln(&fi.elem());

  // Returns a pair with the first being an iterator pointing to the key-value pair and the second a
  // boolean denoting whether a new insertion took place
  auto emplace_ret = _face_to_unc_grad.emplace(&fi, elem_grad);

  mooseAssert(emplace_ret.second, "We should have inserted a new key-value pair");

  VectorValue<ADReal> & unc_face_grad = emplace_ret.first->second;

  const Elem * const neighbor = fi.neighborPtr();

  if (neighbor)
  {
    const VectorValue<ADReal> & neighbor_grad = adGradSln(neighbor);

    const Real g_C = fi.gC();

    // Uncorrected gradient value
    unc_face_grad = g_C * elem_grad + (1. - g_C) * neighbor_grad;
  }
  else // we're on a boundary
  {
    // Do we have any Dirichlet BCs?
    std::vector<FVDirichletBC *> bcs;

    _subproblem.getMooseApp()
        .theWarehouse()
        .query()
        .template condition<AttribSystem>("FVDirichletBC")
        .template condition<AttribThread>(_tid)
        .template condition<AttribBoundaries>(fi.boundaryIDs())
        .template condition<AttribVar>(_var_num)
        .queryInto(bcs);
    mooseAssert(bcs.size() <= 1, "cannot have multiple dirichlet BCs on the same boundary");

    if (bcs.size() == 0)
      // If we don't have a Dirichlet BC, then we should be prescribing some value for the gradient.
      // For now we just prescribe a zero gradient
      unc_face_grad = 0;
    // else
    // We do have a Dirichlet BC. I don't know what we want to do in this case. For now I just use
    // the cell value as the face value
  }

  return unc_face_grad;
}

template <typename OutputType>
const VectorValue<ADReal> &
MooseVariableFV<OutputType>::adGradSln(const FaceInfo & fi) const
{
  auto it = _face_to_grad.find(&fi);

  if (it != _face_to_grad.end())
    return it->second;

  // Returns a pair with the first being an iterator pointing to the key-value pair and the second a
  // boolean denoting whether a new insertion took place
  auto emplace_ret = _face_to_grad.emplace(&fi, uncorrectedAdGradSln(fi));

  mooseAssert(emplace_ret.second, "We should have inserted a new key-value pair");

  VectorValue<ADReal> & face_grad = emplace_ret.first->second;

  const Elem * const neighbor = fi.neighborPtr();

  // perform the correction
  const Point d_CF_vec = neighbor ? fi.neighborCentroid() - fi.elemCentroid()
                                  : (2. * (fi.faceCentroid() - fi.elemCentroid()));
  const Real d_CF = d_CF_vec.norm();

  const Point e_CF = d_CF_vec / d_CF;

  const ADReal elem_value = getElemValue(&fi.elem());

  auto neighbor_value_functor = [&]() {
    if (neighbor)
      return getElemValue(neighbor);
    else
    {
      // If we don't have a neighbor, then we're along a boundary, and we may have a DirichletBC
      std::vector<FVDirichletBC *> bcs;

      _subproblem.getMooseApp()
          .theWarehouse()
          .query()
          .template condition<AttribSystem>("FVDirichletBC")
          .template condition<AttribThread>(_tid)
          .template condition<AttribBoundaries>(fi.boundaryIDs())
          .template condition<AttribVar>(_var_num)
          .queryInto(bcs);
      mooseAssert(bcs.size() <= 1, "cannot have multiple dirichlet BCs on the same boundary");

      bool has_dirichlet_bc = bcs.size() > 0;

      if (has_dirichlet_bc)
      {
        const FVDirichletBC & bc = *bcs[0];

        // Linear interpolation: face_value = (elem_value + neighbor_value) / 2
        return 2. * bc.boundaryValue(fi) - elem_value;
      }
      else
        // No DirichletBC so we'll implicitly apply a zero gradient condition and assume that the
        // face value is equivalent to the element value
        return elem_value;
    }
  };

  face_grad += ((neighbor_value_functor() - elem_value) / d_CF - face_grad * e_CF) * e_CF;

  return face_grad;
}

#endif

template <typename OutputType>
void
MooseVariableFV<OutputType>::residualSetup()
{
  clearCaches();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::jacobianSetup()
{
  clearCaches();
}

template <typename OutputType>
void
MooseVariableFV<OutputType>::clearCaches()
{
#ifdef MOOSE_GLOBAL_AD_INDEXING
  _elem_to_grad.clear();
  _face_to_unc_grad.clear();
  _face_to_grad.clear();
#endif
}

#ifdef MOOSE_GLOBAL_AD_INDEXING
template <typename OutputType>
const ADReal &
MooseVariableFV<OutputType>::adCoeff(const Elem * const elem,
                                     void * context,
                                     ADReal (*fn)(const Elem &, void *)) const
{
  auto it = _elem_to_coeff.find(elem);

  if (it != _elem_to_coeff.end())
    return it->second;

  // Returns a pair with the first being an iterator pointing to the key-value pair and the second a
  // boolean denoting whether a new insertion took place
  auto emplace_ret = _elem_to_coeff.emplace(elem, (*fn)(*elem, context));

  mooseAssert(emplace_ret.second, "We should have inserted a new key-value pair");

  return emplace_ret.first->second;
}
#endif

template class MooseVariableFV<Real>;
// TODO: implement vector fv variable support. This will require some template
// specializations for various member functions in this and the FV variable
// classes. And then you will need to uncomment out the line below:
// template class MooseVariableFV<RealVectorValue>;
