//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "NodeFaceConstraint.h"
#include "MooseEnum.h"

#include "libmesh/coupling_matrix.h"

#include <vector>
#include <unordered_map>

class RANFSNormalMechanicalContact;
class PenetrationInfo;

namespace libMesh
{
template <typename>
class NumericVector;
}

template <>
InputParameters validParams<RANFSNormalMechanicalContact>();

class RANFSNormalMechanicalContact : public NodeFaceConstraint
{
public:
  RANFSNormalMechanicalContact(const InputParameters & parameters);

  bool shouldApply() override;
  void residualSetup() override;
  void timestepSetup() override;
  void initialSetup() override;
  bool overwriteSlaveResidual() override;
  void computeSlaveValue(NumericVector<Number> & solution) override;

protected:
  Real computeQpSlaveValue() override;

  Real computeQpResidual(Moose::ConstraintType type) override;
  Real computeQpJacobian(Moose::ConstraintJacobianType type) override;
  Real computeQpOffDiagJacobian(Moose::ConstraintJacobianType type, unsigned int jvar) override;

  const MooseEnum _component;
  const unsigned int _mesh_dimension;
  NumericVector<Number> & _residual_copy;

  unsigned int _largest_component;
  std::vector<unsigned int> _vars;
  std::vector<MooseVariable *> _var_objects;
  std::unordered_map<dof_id_type, Real> _node_to_contact_lm;
  std::unordered_map<dof_id_type, Real> _node_to_tied_lm;
  std::unordered_map<dof_id_type, std::vector<const Elem *>> _node_to_master_elem_sequence;
  Real _contact_lm;
  Real _tied_lm;
  PenetrationInfo * _pinfo;
  std::unordered_map<dof_id_type, const Node *> _ping_pong_slave_node_to_master_node;
  Real _distance;
  bool _tie_nodes;
  unsigned int _master_index;
  RealVectorValue _res_vec;
  const Node * _nearest_node;

  /// Vector size corresponds to the number of displacement variables. Each index in the vector
  /// holds an unordered map. The map container maps a non-zero global column index to its
  /// corresponding entry in the Jacobian. This container is very useful in conjunction with
  /// NodeFaceConstraint::_connected_dof_indices which is a vector whose entries correspond to the
  /// global degrees of freedom of the current jvar
  std::vector<std::unordered_map<dof_id_type, Number>> _dof_number_to_value;

  CouplingMatrix _disp_coupling;
  const bool _ping_pong_protection;

  std::unordered_map<unsigned int, unsigned int> _jvar_to_comp;
};
