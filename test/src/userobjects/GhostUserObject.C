//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GhostUserObject.h"
#include "MooseMesh.h"

// invalid_processor_id
#include "libmesh/dof_object.h"

registerMooseObject("MooseTestApp", GhostUserObject);

template <>
InputParameters
validParams<GhostUserObject>()
{
  InputParameters params = validParams<ElementUserObject>();
  params.addParam<unsigned int>(
      "rank",
      DofObject::invalid_processor_id,
      "The rank for which the ghosted elements are recorded (Default: ALL)");

  //  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;

  params.registerRelationshipManagers("ElementSideNeighborLayers");
  params.addRequiredParam<unsigned short>("element_side_neighbor_layers",
                                          "Number of layers to ghost");

  params.addClassDescription("User object to calculate ghosted elements on a single processor or "
                             "the union across all processors.");

  params.addRequiredCoupledVar(
      "some_variable", "The gradient of this variable will be used as the velocity vector.");

  return params;
}

GhostUserObject::GhostUserObject(const InputParameters & parameters)
  : ElementUserObject(parameters),
    _rank(getParam<unsigned int>("rank")),
    _some_variable(coupledValue("some_variable"))
{
}

void
GhostUserObject::initialize()
{
  _ghost_data.clear();
}

void
GhostUserObject::execute()
{
  const auto & mesh = _subproblem.mesh().getMesh();

  auto my_processor_id = processor_id();

  for (const auto & elem : mesh.active_element_ptr_range())
  {
    _console << "Elem ID: " << elem->id() << std::flush;

    _fe_problem.prepare(elem, _tid);
    _fe_problem.reinitElem(elem, _tid);

    if (my_processor_id == 0 && elem->processor_id() != my_processor_id)
    {
      _console << ": " << _some_variable[0] << std::endl;
      _ghost_data[elem->id()] = _some_variable[0];
    }

    //
    //
    //  if (_rank == DofObject::invalid_processor_id || my_processor_id == _rank)
    //  {
    //    const auto & mesh = _subproblem.mesh().getMesh();
    //
    //    for (const auto & elem : mesh.active_element_ptr_range())
    //      if (elem->processor_id() != my_processor_id)
    //        _ghost_data.emplace(elem->id());
    //  }
  }
}

void
GhostUserObject::threadJoin(const UserObject &)
{
}

void
GhostUserObject::finalize()
{
  _communicator.set_union(_ghost_data);
}

Real
GhostUserObject::getElementalValue(dof_id_type element_id) const
{
  auto pos = _ghost_data.find(element_id);

  if (pos != _ghost_data.end())
    return pos->second;
  else
    return 0;
}
