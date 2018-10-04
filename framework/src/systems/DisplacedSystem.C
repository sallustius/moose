//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DisplacedSystem.h"
#include "DisplacedProblem.h"
#include "MaxVarNDofsPerElem.h"

DisplacedSystem::DisplacedSystem(DisplacedProblem & problem,
                                 SystemBase & undisplaced_system,
                                 const std::string & name,
                                 Moose::VarKindType var_kind)
  : SystemBase(problem, name, var_kind),
    _fe_problem(problem.feProblem()),
    _undisplaced_system(undisplaced_system),
    _sys(problem.es().add_system<TransientExplicitSystem>(name))
{
}

DisplacedSystem::~DisplacedSystem() {}

void
DisplacedSystem::init()
{
  Moose::perf_log.push("maxVarNDofsPerElem()", "Setup");
  MaxVarNDofsPerElem mvndpe(_fe_problem, _undisplaced_system);
  Threads::parallel_reduce(*_fe_problem.mesh().getActiveLocalElementRange(), mvndpe);
  _max_var_n_dofs_per_elem = mvndpe.max();
  _communicator.max(_max_var_n_dofs_per_elem);
  Moose::perf_log.pop("maxVarNDofsPerElem()", "Setup");
}

NumericVector<Number> &
DisplacedSystem::getVector(const std::string & name)
{
  if (_sys.have_vector(name))
    return _sys.get_vector(name);
  else
    return _undisplaced_system.getVector(name);
}
