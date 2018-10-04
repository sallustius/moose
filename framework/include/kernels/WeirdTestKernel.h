//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef WEIRDTESTKERNEL_H
#define WEIRDTESTKERNEL_H

#include "Kernel.h"

class WeirdTestKernel;

template <>
InputParameters validParams<WeirdTestKernel>();

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class WeirdTestKernel : public Kernel
{
public:
  WeirdTestKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;

  const unsigned int _disp_x_id;
  const unsigned int _disp_y_id;
  const unsigned int _disp_z_id;
  const unsigned int _dim;
  const MooseArray<std::vector<NumberArray<AD_MAX_DOFS_PER_ELEM, Real>>> & _dphidx_derivatives;
  const MooseArray<ADRealGradient> & _ad_grad_u;
};

#endif /* WEIRDTESTKERNEL_H */
