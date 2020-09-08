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
#include "FVUtils.h"

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

template <typename OutputType>
std::pair<bool, const FVDirichletBC *>
MooseVariableFV<OutputType>::getDirichletBC(const FaceInfo & fi) const
{
  std::vector<FVDirichletBC *> bcs;

  _subproblem.getMooseApp()
      .theWarehouse()
      .query()
      .template condition<AttribSystem>("FVDirichletBC")
      .template condition<AttribThread>(_tid)
      .template condition<AttribBoundaries>(fi.boundaryIDs())
      .template condition<AttribVar>(_var_num)
      .template condition<AttribSysNum>(_sys.number())
      .queryInto(bcs);
  mooseAssert(bcs.size() <= 1, "cannot have multiple dirichlet BCs on the same boundary");

  bool has_dirichlet_bc = bcs.size() > 0;

  if (has_dirichlet_bc)
  {
    mooseAssert(bcs.size() == 1,
                "There should not be multiple Dirichlet boundary conditions for a given variable "
                "on a given FaceInfo");

    mooseAssert(bcs[0], "The FVDirichletBC is null!");

    return std::make_pair(true, bcs[0]);
  }
  else
    return std::make_pair(false, nullptr);
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
ADReal
MooseVariableFV<OutputType>::getNeighborValue(const Elem * const neighbor,
                                              const FaceInfo & fi,
                                              const ADReal & elem_value) const
{
  if (neighbor)
    return getElemValue(neighbor);
  else
  {
    // If we don't have a neighbor, then we're along a boundary, and we may have a DirichletBC
    const auto & pr = getDirichletBC(fi);

    if (pr.first)
    {
      mooseAssert(pr.second, "The FVDirichletBC is null!");

      const FVDirichletBC & bc = *pr.second;

      // Linear interpolation: face_value = (elem_value + neighbor_value) / 2
      return 2. * bc.boundaryValue(fi) - elem_value;
    }
    else
      // No DirichletBC so we'll implicitly apply a zero gradient condition and assume that the
      // face value is equivalent to the element value
      return elem_value;
  }
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

  auto action_functor = [&grad, &volume_set, &volume, &elem_value, this](
                            const Elem & functor_elem,
                            const Elem * const neighbor,
                            const FaceInfo * const fi,
                            const Point & surface_vector,
                            Real coord,
                            const bool elem_has_info) {
    auto face_value_functor = [&neighbor, &fi, &elem_value, this]() {
      if (neighbor)
      {
        ADReal neighbor_value = getElemValue(neighbor);

        const Real g_C = fi->gC();

        return g_C * elem_value + (1. - g_C) * neighbor_value;
      }
      else
      {
        // If we don't have a neighbor, then we're along a boundary, and we may have a DirichletBC
        const auto & pr = getDirichletBC(*fi);

        if (pr.first)
        {
          const FVDirichletBC & bc = *pr.second;

          return ADReal(bc.boundaryValue(*fi));
        }
        else
        {
          // No DirichletBC so we'll implicitly apply a zero gradient condition and assume that the
          // face value is equivalent to the element value
          return elem_value;
        }
      }
    };

    grad += face_value_functor() * surface_vector;

    if (!volume_set)
    {
      if (elem_has_info)
      {
        coordTransformFactor(_subproblem, functor_elem.subdomain_id(), fi->elemCentroid(), coord);
        volume = fi->elemVolume() * coord;
      }
      else
      {
        coordTransformFactor(_subproblem, neighbor->subdomain_id(), fi->neighborCentroid(), coord);
        volume = fi->neighborVolume() * coord;
      }

      volume_set = true;
    }
  };

  Moose::loopOverElemFaceInfo(*elem, _mesh, _subproblem, action_functor);

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

  // If we have a neighbor then we interpolate between the two to the face. If we do not, then we
  // check for a Dirichlet BC. If we have a Dirichlet BC, then we will apply a zero Hessian
  // assumption. If we do not, then we know we are applying a zero gradient assumption elsehwere in
  // our calculations, so we should be consistent and apply a zero gradient assumption here as well
  if (neighbor)
  {
    const VectorValue<ADReal> & neighbor_grad = adGradSln(neighbor);

    const Real g_C = fi.gC();

    // Uncorrected gradient value
    unc_face_grad = g_C * elem_grad + (1. - g_C) * neighbor_grad;
  }
  else
  {
    const auto & pr = getDirichletBC(fi);

    if (!pr.first)
      unc_face_grad = 0;
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

  face_grad +=
      ((getNeighborValue(neighbor, fi, elem_value) - elem_value) / d_CF - face_grad * e_CF) * e_CF;

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
                                     ADReal (*fn)(const Elem * const, void *)) const
{
  auto it = _elem_to_coeff.find(elem);

  if (it != _elem_to_coeff.end())
    return it->second;

  // Returns a pair with the first being an iterator pointing to the key-value pair and the second a
  // boolean denoting whether a new insertion took place
  auto emplace_ret = _elem_to_coeff.emplace(elem, (*fn)(elem, context));

  mooseAssert(emplace_ret.second, "We should have inserted a new key-value pair");

  return emplace_ret.first->second;
}
#endif

template class MooseVariableFV<Real>;
// TODO: implement vector fv variable support. This will require some template
// specializations for various member functions in this and the FV variable
// classes. And then you will need to uncomment out the line below:
// template class MooseVariableFV<RealVectorValue>;
