//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INSADMOMENTUMLAPLACEFORM_H
#define INSADMOMENTUMLAPLACEFORM_H

#include "INSADMomentumBase.h"

// Forward Declarations
template <ComputeStage>
class INSADMomentumLaplaceForm;

declareADValidParams(INSADMomentumLaplaceForm);

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "Laplacian" form of the governing equations.
 */
template <ComputeStage compute_stage>
class INSADMomentumLaplaceForm : public INSADMomentumBase<compute_stage>
{
public:
  INSADMomentumLaplaceForm(const InputParameters & parameters);

  virtual ~INSADMomentumLaplaceForm() {}

protected:
  virtual ADResidual computeQpResidualViscousPart() override;

  usingINSMomentumBaseMembers;
};

#endif
