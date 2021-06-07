//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADReal.h"
#include <memory>

class MooseEnum;

namespace Moose
{
namespace FV
{
enum class LimiterType : int
{
  VanLeer = 0,
  Upwind,
  CentralDifference,
  MinMod,
  SOU,
  QUICK
};
extern const MooseEnum moose_limiter_type;

/**
 * Base class for defining slope limiters for finite volume or potentially reconstructed
 * Discontinuous-Galerkin applications
 */
class Limiter
{
public:
  /**
   * Defines the slope limiter function $\beta(r_f)$ where $r_f$ represents the ratio of upstream to
   * downstream gradients
   */
  virtual ADReal operator()(const ADReal & r_f) const = 0;

  Limiter() = default;

  virtual ~Limiter() = default;

  static std::unique_ptr<Limiter> build(LimiterType limiter);
};

}
}
