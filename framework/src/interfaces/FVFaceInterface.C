//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVFaceInterface.h"
#include "MooseObject.h"
#include "InputParameters.h"

FVFaceInterface::FVFaceInterface(const MooseObject * const moose_object)
  : _fvfi_subproblem(
        *moose_object->parameters().getCheckedPointerParam<SubProblem *>("_subproblem"))
{
}
