//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Limiter.h"

namespace Moose
{
namespace FV
{
/**
 * Implements a limiter which reproduces the upwind scheme, defined by
 * $\beta(r_f) = 0$
 */
class UpwindLimiter : public Limiter
{
public:
  ADReal operator()(const ADReal & /*r_f*/) const override final { return 0; }

  UpwindLimiter() = default;
};
}
}
