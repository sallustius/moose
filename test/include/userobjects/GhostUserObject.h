//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef GHOSTUSEROBJECT_H
#define GHOSTUSEROBJECT_H

#include "GeneralUserObject.h"
#include "Coupleable.h"
#include "MooseVariableDependencyInterface.h"

// Forward Declarations
class GhostUserObject;

template <>
InputParameters validParams<GhostUserObject>();

/**
 * User object to calculate ghosted elements on a single processor or the union across all
 * processors.
 */
class GhostUserObject : public GeneralUserObject,
                        public Coupleable,
                        public MooseVariableDependencyInterface
{
public:
  GhostUserObject(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject &) override;

  Real getElementalValue(dof_id_type element_id) const;

protected:
  void buildRange();

  std::vector<MooseVariableFEBase *> _fe_vars;

  std::map<dof_id_type, Real> _ghost_data;

  dof_id_type _rank;

  const VariableValue & _some_variable;

  std::unique_ptr<StoredRange<std::set<const Elem *>::iterator, const Elem *>> _elem_range;
};

#endif // GHOSTUSEROBJECT_H
