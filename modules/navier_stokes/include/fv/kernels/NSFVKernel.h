//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVMatAdvection.h"

#ifdef MOOSE_GLOBAL_AD_INDEXING

template <typename>
class MooseVariableFV;

namespace libMesh
{
class DofMap;
template <typename>
class NumericVector;
template <typename>
class VectorValue;
}

class NSFVKernel : public FVMatAdvection
{
public:
  static InputParameters validParams();
  NSFVKernel(const InputParameters & params);

  ADReal coeffCalculator(const Elem * elem);

protected:
  /**
   * interpolation overload for the velocity
   */
  void interpolate(InterpMethod m,
                   ADRealVectorValue & interp_v,
                   const ADRealVectorValue & elem_v,
                   const ADRealVectorValue & neighbor_v);

  ADReal computeQpResidual() override;

private:
  /// pressure variable
  const MooseVariableFV<Real> * const _p_var;
  /// x-velocity
  const MooseVariableFV<Real> * const _u_var;
  /// y-velocity
  const MooseVariableFV<Real> * const _v_var;
  /// z-velocity
  const MooseVariableFV<Real> * const _w_var;

  /// The viscosity
  const Real _mu;
  /// The density
  const Real _rho;

  /// The interpolation method to use for the velocity
  InterpMethod _velocity_interp_method;
};

#endif
