// App includes
#include "AugmentSparsityOnInterface.h"
#include "Executioner.h"
#include "FEProblemBase.h"
#include "MooseApp.h"

// libMesh includes
#include "libmesh/elem.h"

registerMooseObject("MooseApp", AugmentSparsityOnInterface);

using namespace libMesh;

template <>
InputParameters
validParams<AugmentSparsityOnInterface>()
{
  InputParameters params = validParams<AlgebraicRelationshipManager>();
  params.addRequiredParam<SubdomainID>("master_subdomain_id", "The id of the master subdomain.");
  params.addRequiredParam<SubdomainID>("slave_subdomain_id", "The id of the slave subdomain.");
  return params;
}

AugmentSparsityOnInterface::AugmentSparsityOnInterface(const InputParameters & params)
  : AlgebraicRelationshipManager(params),
    _master_subdomain_id(getParam<SubdomainID>("master_subdomain_id")),
    _slave_subdomain_id(getParam<SubdomainID>("slave_subdomain_id"))
{
  _rm_type = Moose::RelationshipManagerType::ALGEBRAIC;
}

void
AugmentSparsityOnInterface::mesh_reinit()
{
  auto n_active_elem = _mesh.getMesh().n_active_elem();

  for (unsigned int i = 0; i < n_active_elem; ++i)
  {
    Elem * elem = _mesh.getMesh().elem_ptr(i);
    if (elem->subdomain_id() == _slave_subdomain_id || elem->subdomain_id() == _master_subdomain_id)
    {
      _interface_elems.insert(elem);
      _interface_elems.insert(elem->interior_parent());
    }
  }

  // for (MeshBase::const_element_iterator el = _mesh.getMesh().active_elements_begin(),
  //                                       end_el = _mesh.getMesh().active_elements_end();
  //      el != end_el;
  //      ++el)
  // {
  //   const Elem * elem = *el;
  //   if (elem->subdomain_id() == _slave_subdomain_id || elem->subdomain_id() ==
  //   _master_subdomain_id)
  //   {
  //     _interface_elems.insert(elem);
  //     _interface_elems.insert(elem->interior_parent());
  //   }
  // }
}

void
AugmentSparsityOnInterface::attachRelationshipManagersInternal(
    Moose::RelationshipManagerType rm_type)
{
  if (rm_type == Moose::RelationshipManagerType::GEOMETRIC)
    attachGeometricFunctorHelper(*this);
  else
    attachAlgebraicFunctorHelper(*this);
}

std::string
AugmentSparsityOnInterface::getInfo() const
{
  std::ostringstream oss;
  oss << "AugmentSparsityOnInterface";
  return oss.str();
}

void
AugmentSparsityOnInterface::operator()(const MeshBase::const_element_iterator & range_begin,
                                       const MeshBase::const_element_iterator & range_end,
                                       processor_id_type p,
                                       map_type & coupled_elements)
{
  const CouplingMatrix * const null_mat = libmesh_nullptr;

  for (const auto & elem : as_range(range_begin, range_end))
    if (elem->processor_id() != p && _interface_elems.find(elem) != _interface_elems.end())
    {
      coupled_elements.insert(std::make_pair(elem, null_mat));
      std::cout << elem->centroid() << std::endl;
    }
}
