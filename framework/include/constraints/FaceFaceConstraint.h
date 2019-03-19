#ifndef FACEFACECONSTRAINT_H
#define FACEFACECONSTRAINT_H

#include "MortarConstraint.h"

template <ComputeStage>
class FaceFaceConstraint;

declareADValidParams(FaceFaceConstraint);

/**
 * This is a deprecated object!  Use MortarConstraint instead!
 */
template <ComputeStage compute_stage>
class FaceFaceConstraint : public MortarConstraint<compute_stage>
{
public:
  FaceFaceConstraint(const InputParameters & params) : MortarConstraint<compute_stage>(params)
  {
    mooseDeprecated("FaceFaceConstraint is deprecated!  Use MortarConstraint instead!");
  }
};

#endif
