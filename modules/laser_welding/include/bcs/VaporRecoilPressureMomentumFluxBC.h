//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef VAPORRECOILPRESSUREMOMENTUMFLUXBC_H
#define VAPORRECOILPRESSUREMOMENTUMFLUXBC_H

#include "ADIntegratedBC.h"

template <ComputeStage compute_stage>
class VaporRecoilPressureMomentumFluxBC;

declareADValidParams(VaporRecoilPressureMomentumFluxBC);

template <ComputeStage compute_stage>
class VaporRecoilPressureMomentumFluxBC : public ADIntegratedBC<compute_stage>
{
public:
  VaporRecoilPressureMomentumFluxBC(const InputParameters & parameters);

protected:
  virtual ADResidual computeQpResidual() override;

  usingIntegratedBCMembers;

  const Real _ap0;
  const Real _ap1;
  const Real _ap2;
  const Real _ap3;
  const Real _bp0;
  const Real _bp1;
  const Real _bp2;
  const Real _bp3;
  const Real _Tb;
  const Real _Tbound1;
  const Real _Tbound2;

  const ADVariableValue & _temperature;
  const unsigned _component;
};

#endif /* VAPORRECOILPRESSUREMOMENTUMFLUXBC_H */
