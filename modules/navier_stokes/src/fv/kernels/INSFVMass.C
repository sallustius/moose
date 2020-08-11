//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMass.h"

#ifdef MOOSE_GLOBAL_AD_INDEXING

#include "NS.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "MooseMesh.h"

#include "libmesh/dof_map.h"
#include "libmesh/mesh_base.h"

registerMooseObject("NavierStokesApp", INSFVMass);

InputParameters
INSFVMass::validParams()
{
  InputParameters params = NSFVKernel::validParams();
  params.addClassDescription("Enforces the divergence free condition of the velocity field.");
  params.set<MaterialPropertyName>("advected_quantity") = std::to_string(1);
  params.set<MaterialPropertyName>("vel") = NS::velocity;
  params.addRequiredParam<bool>(
      "constrain_pressure",
      "Whether to perform a pressure pin. This is required at least for problems where Dirichlet "
      "values are prescribed for all normal velocity components for all boundaries.");
  return params;
}

INSFVMass::INSFVMass(const InputParameters & params)
  : NSFVKernel(params), _constrain_pressure(getParam<bool>("constrain_pressure"))
{
}

void
INSFVMass::initialSetup()
{
  if (_constrain_pressure)
  {
    DofMap & dof_map = _subproblem.systemBaseNonlinear().dofMap();

    // No default constructors for these unfortunately
    auto begin_elem_it = _mesh.getMesh().active_local_elements_begin();
    auto end_elem_it = _mesh.getMesh().active_local_elements_end();

    if (blockRestricted())
    {
      mooseAssert(blockIDs().begin() != blockIDs().end(), "blockIDs shouldn't be empty");
      begin_elem_it = _mesh.getMesh().active_local_subdomain_elements_begin(*blockIDs().begin());
      end_elem_it = _mesh.getMesh().active_local_subdomain_elements_end(*blockIDs().begin());
    }

    std::vector<bool> global_found_element(_communicator.size(), false);
    global_found_element[_communicator.rank()] = (begin_elem_it != end_elem_it);

    _communicator.max(global_found_element);

    auto it = std::find(global_found_element.begin(), global_found_element.end(), true);

    mooseAssert(it != global_found_element.end(),
                "We had to have found an element matching this kernel somewhere!");

    bool apply_my_constraint =
        std::distance(global_found_element.begin(), it) == _communicator.rank();

    if (!apply_my_constraint)
      return;

    std::vector<dof_id_type> pressure_dof_indices;
    dof_map.dof_indices(*begin_elem_it, pressure_dof_indices, _var.number());

    mooseAssert(pressure_dof_indices.size(), "We should have gotten some dof indices");

    const auto dof_to_constrain = pressure_dof_indices.front();

    dof_map.add_constraint_row(dof_to_constrain, DofConstraintRow{}, 0, true);
  }
}

#endif
