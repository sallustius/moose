//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressDivergenceTensors.h"

// MOOSE includes
#include "ElasticityTensorTools.h"
#include "Material.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

registerMooseObject("TensorMechanicsApp", StressDivergenceTensors);

template <>
InputParameters
validParams<StressDivergenceTensors>()
{
  InputParameters params = validParams<ADKernel>();
  params.addClassDescription("Stress divergence kernel for the Cartesian coordinate system");
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in. (0 for x, "
                                        "1 for y, 2 for z)");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  params.addCoupledVar("temperature",
                       "The name of the temperature variable used in the "
                       "ComputeThermalExpansionEigenstrain.  (Not required for "
                       "simulations without temperature coupling.)");
  params.addParam<std::string>(
      "thermal_eigenstrain_name",
      "thermal_eigenstrain",
      "The eigenstrain_name used in the ComputeThermalExpansionEigenstrain.");
  params.addCoupledVar("out_of_plane_strain",
                       "The name of the out_of_plane_strain variable used in the "
                       "WeakPlaneStress kernel. Required only if want to provide off-diagonal "
                       "Jacobian in plane stress analysis using weak formulation.");
  MooseEnum out_of_plane_direction("x y z", "z");
  params.addParam<MooseEnum>(
      "out_of_plane_direction",
      out_of_plane_direction,
      "The direction of the out_of_plane_strain variable used in the WeakPlaneStress kernel.");
  params.addParam<std::string>("base_name", "Material property base name");
  params.set<bool>("use_displaced_mesh") = false;
  params.addParam<bool>(
      "use_finite_deform_jacobian", false, "Jacobian for corotational finite strain");
  params.addParam<bool>("volumetric_locking_correction",
                        false,
                        "Set to false to turn off volumetric locking correction");
  return params;
}

StressDivergenceTensors::StressDivergenceTensors(const InputParameters & parameters)
  : ADKernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _stress(getADMaterialPropertyByName<RankTwoTensor>(_base_name + "stress")),
    _Jacobian_mult(getADMaterialPropertyByName<RankFourTensor>(_base_name + "Jacobian_mult")),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _temp_coupled(isCoupled("temperature")),
    _temp_var(_temp_coupled ? coupled("temperature") : 0),
    _out_of_plane_strain_coupled(isCoupled("out_of_plane_strain")),
    _out_of_plane_strain_var(_out_of_plane_strain_coupled ? coupled("out_of_plane_strain") : 0),
    _out_of_plane_direction(getParam<MooseEnum>("out_of_plane_direction")),
    _avg_grad_test(_test.size(), std::vector<Real>(3, 0.0)),
    _avg_grad_phi(_phi.size(), std::vector<Real>(3, 0.0)),
    _volumetric_locking_correction(getParam<bool>("volumetric_locking_correction"))
{
  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

  // Checking for consistency between mesh size and length of the provided displacements vector
  if (_out_of_plane_direction != 2 && _ndisp != 3)
    mooseError("For 2D simulations where the out-of-plane direction is x or y coordinate "
               "directions the number of supplied displacements must be three.");
  else if (_out_of_plane_direction == 2 && _ndisp != _mesh.dimension())
    mooseError("The number of displacement variables supplied must match the mesh dimension");

  // Error if volumetic locking correction is turned on for 1D problems
  if (_ndisp == 1 && _volumetric_locking_correction)
    mooseError("Volumetric locking correction should be set to false for 1-D problems.");
}

void
StressDivergenceTensors::initialSetup()
{
  if (getBlockCoordSystem() != Moose::COORD_XYZ)
    mooseError(
        "The coordinate system in the Problem block must be set to XYZ for cartesian geometries.");
}

void
StressDivergenceTensors::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  _local_re.resize(re.size());
  _local_re.zero();

  if (_volumetric_locking_correction)
    computeAverageGradientTest();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual().value();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (const auto & var : _save_in)
      var->sys().solution().add_vector(_local_re, var->dofIndices());
  }
}

ADReal
StressDivergenceTensors::computeQpResidual()
{
  Real residual = _stress[_qp].row(_component) * _grad_test[_i][_qp];
  // volumetric locking correction
  if (_volumetric_locking_correction)
    residual += _stress[_qp].tr() / 3.0 *
                (_avg_grad_test[_i][_component] - _grad_test[_i][_qp](_component));

  return residual;
}

void
StressDivergenceTensors::computeJacobian()
{
  if (_volumetric_locking_correction)
    computeAverageGradientTest();
  ADKernel::computeJacobian();
}

void
StressDivergenceTensors::computeOffDiagJacobian(MooseVariableFEBase & jvar)
{
  if (_volumetric_locking_correction)
    computeAverageGradientTest();
  ADKernel::computeOffDiagJacobian(jvar);
}

void
StressDivergenceTensors::computeAverageGradientTest()
{
  // Calculate volume averaged value of shape function derivative
  _avg_grad_test.resize(_test.size());
  for (_i = 0; _i < _test.size(); ++_i)
  {
    _avg_grad_test[_i].resize(3);
    _avg_grad_test[_i][_component] = 0.0;
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      _avg_grad_test[_i][_component] += _grad_test[_i][_qp](_component) * _JxW[_qp] * _coord[_qp];

    _avg_grad_test[_i][_component] /= _current_elem_volume;
  }
}
