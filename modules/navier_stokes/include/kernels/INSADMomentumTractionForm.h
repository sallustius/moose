//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INSADMOMENTUMTRACTIONFORM_H
#define INSADMOMENTUMTRACTIONFORM_H

#include "INSADMomentumBase.h"

// Forward Declarations
template <ComputeStage>
class INSADMomentumTractionForm;

declareADValidParams(INSADMomentumTractionForm);

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "traction" form of the governing equations.
 */
template <ComputeStage compute_stage>
class INSADMomentumTractionForm : public INSADMomentumBase<compute_stage>
{
public:
  INSADMomentumTractionForm(const InputParameters & parameters);

  virtual ~INSADMomentumTractionForm() {}

protected:
  virtual ADResidual computeQpResidualViscousPart() override;

  usingINSMomentumBaseMembers;
};

#endif
