//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADMass.h"
#include "Function.h"

registerADMooseObject("NavierStokesApp", INSADMass);

defineADValidParams(
    INSADMass,
    INSADBase,
    params.addClassDescription("This class computes the mass equation residual and Jacobian "
                               "contributions for the incompressible Navier-Stokes momentum "
                               "equation.");
    params.addParam<bool>("pspg",
                          false,
                          "Whether to perform PSPG stabilization of the mass equation");
    params.addParam<FunctionName>("x_vel_forcing_func", 0, "The x-velocity mms forcing function.");
    params.addParam<FunctionName>("y_vel_forcing_func", 0, "The y-velocity mms forcing function.");
    params.addParam<FunctionName>("z_vel_forcing_func",
                                  0,
                                  "The z-velocity mms forcing function."););

template <ComputeStage compute_stage>
INSADMass<compute_stage>::INSADMass(const InputParameters & parameters)
  : INSADBase<compute_stage>(parameters),
    _pspg(adGetParam<bool>("pspg")),
    _x_ffn(this->getFunction("x_vel_forcing_func")),
    _y_ffn(this->getFunction("y_vel_forcing_func")),
    _z_ffn(this->getFunction("z_vel_forcing_func"))

{
}

template <ComputeStage compute_stage>
ADResidual
INSADMass<compute_stage>::computeQpResidual()
{
  // (div u) * q
  // Note: we (arbitrarily) multiply this term by -1 so that it matches the -p(div v)
  // term in the momentum equation.  Not sure if that is really important?
  auto r = -(_grad_u_vel[_qp](0) + _grad_v_vel[_qp](1) + _grad_w_vel[_qp](2)) * _test[_i][_qp];

  if (_pspg)
    r += computeQpPGResidual();

  return r;
}

template <ComputeStage compute_stage>
ADResidual
INSADMass<compute_stage>::computeQpPGResidual()
{
  const auto & viscous_term =
      _laplace ? this->strongViscousTermLaplace() : this->strongViscousTermTraction();
  const auto & transient_term =
      _transient_term ? this->timeDerivativeTerm() : INSVectorValue<compute_stage>(0, 0, 0);
  const auto & convective_term =
      _convective_term ? this->convectiveTerm() : INSVectorValue<compute_stage>(0, 0, 0);
  auto r = -1. / _rho[_qp] * this->tau() * _grad_test[_i][_qp] *
           (this->strongPressureTerm() + this->gravityTerm() + viscous_term + convective_term +
            transient_term -
            RealVectorValue(_x_ffn.value(_t, _q_point[_qp]),
                            _y_ffn.value(_t, _q_point[_qp]),
                            _z_ffn.value(_t, _q_point[_qp])));

  return r;
}

template class INSADMass<RESIDUAL>;
template class INSADMass<JACOBIAN>;
