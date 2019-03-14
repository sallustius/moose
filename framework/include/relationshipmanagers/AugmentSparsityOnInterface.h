#ifndef AUGMENT_SPARSITY_ON_INTERFACE_H
#define AUGMENT_SPARSITY_ON_INTERFACE_H

// App includes
#include "AutomaticMortarGeneration.h"
#include "AlgebraicRelationshipManager.h"

// libMesh includes
#include "libmesh/mesh_base.h"

using libMesh::boundary_id_type;
using libMesh::Elem;
using libMesh::GhostingFunctor;
using libMesh::MeshBase;
using libMesh::processor_id_type;

class AugmentSparsityOnInterface;

template <>
InputParameters validParams<AugmentSparsityOnInterface>();

class AugmentSparsityOnInterface : public AlgebraicRelationshipManager
{
public:
  AugmentSparsityOnInterface(const InputParameters &);

  /**
   * This function must be overriden by application codes to add
   * required elements from (range_begin, range_end) to the
   * coupled_elements map.
   */
  virtual void operator()(const MeshBase::const_element_iterator & range_begin,
                          const MeshBase::const_element_iterator & range_end,
                          processor_id_type p,
                          map_type & coupled_elements) override;

  /**
   * According to the base class docs, "We call mesh_reinit() whenever
   * the relevant Mesh has changed, but before remote elements on a
   * distributed mesh are deleted."
   */
  virtual void mesh_reinit() override;

  /**
   * Update the cached _lower_to_upper map whenever our Mesh has been
   * redistributed.  We'll be lazy and just recalculate from scratch.
   */
  virtual void redistribute() override { this->mesh_reinit(); }

  void attachRelationshipManagersInternal(Moose::RelationshipManagerType rm_type) override;

  std::string getInfo() const override;

protected:
  SubdomainID _master_subdomain_id;
  SubdomainID _slave_subdomain_id;

  std::set<const Elem *> _interface_elems;
};

#endif
