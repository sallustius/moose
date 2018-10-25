//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef INSADTEMPERATURE_H
#define INSADTEMPERATURE_H

#include "ADKernel.h"

// Forward Declarations
template <ComputeStage>
class INSADTemperature;

declareADValidParams(INSADTemperature);

/**
 * This class computes the residual and Jacobian contributions for the
 * incompressible Navier-Stokes temperature (energy) equation.
 */
template <ComputeStage compute_stage>
class INSADTemperature : public ADKernel<compute_stage>
{
public:
  INSADTemperature(const InputParameters & parameters);

  virtual ~INSADTemperature() {}

protected:
  virtual ADResidual computeQpResidual() override;

  // Coupled variables
  const ADVariableValue & _u_vel;
  const ADVariableValue & _v_vel;
  const ADVariableValue & _w_vel;

  // Required parameters
  const ADMaterialProperty(Real) & _rho;
  const ADMaterialProperty(Real) & _k;
  const ADMaterialProperty(Real) & _cp;

  usingKernelMembers;
};

#endif // INSADTEMPERATURE_H
