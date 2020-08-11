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

  ADReal coeffCalculator(const Elem & elem);

protected:
  ADRealVectorValue interpolateVelocity();

  ADReal computeQpResidual() override;

private:
  const MooseVariableFV<Real> * const _p_var;
  const MooseVariableFV<Real> * const _u_var;
  const MooseVariableFV<Real> * const _v_var;
  const MooseVariableFV<Real> * const _w_var;

  const Real _mu;
  const Real _rho;
};

#endif
