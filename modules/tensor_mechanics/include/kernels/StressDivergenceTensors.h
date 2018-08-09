//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef STRESSDIVERGENCETENSORS_H
#define STRESSDIVERGENCETENSORS_H

#include "ADKernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"

// Forward Declarations
class StressDivergenceTensors;
class RankTwoTensor;
class RankFourTensor;

template <>
InputParameters validParams<StressDivergenceTensors>();

/**
 * StressDivergenceTensors mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class StressDivergenceTensors : public ADKernel
{
public:
  StressDivergenceTensors(const InputParameters & parameters);

  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
  using ADKernel::computeOffDiagJacobian;

protected:
  virtual void initialSetup() override;

  virtual void computeResidual() override;
  virtual ADReal computeQpResidual() override;

  virtual void computeAverageGradientTest();
  virtual void computeAverageGradientPhi();

  std::string _base_name;

  const ADMaterialProperty<RankTwoTensor> & _stress;
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;

  const unsigned int _component;

  /// Coupled displacement variables
  unsigned int _ndisp;
  std::vector<unsigned int> _disp_var;

  const bool _temp_coupled;
  const unsigned int _temp_var;

  const bool _out_of_plane_strain_coupled;
  const unsigned int _out_of_plane_strain_var;
  const unsigned int _out_of_plane_direction;

  /// Gradient of test function averaged over the element. Used in volumetric locking correction calculation.
  std::vector<std::vector<Real>> _avg_grad_test;

  /// Gradient of phi function averaged over the element. Used in volumetric locking correction calculation.
  std::vector<std::vector<Real>> _avg_grad_phi;

  /// Flag for volumetric locking correction
  bool _volumetric_locking_correction;
};

#endif // STRESSDIVERGENCETENSORS_H
