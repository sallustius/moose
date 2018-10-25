//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSADMomentumBase.h"
#include "Function.h"

defineADValidParams(
    INSADMomentumBase,
    INSADBase,
    params.addRequiredParam<unsigned>("component",
                                      "The velocity component that this is applied to.");
    params.addParam<bool>("integrate_p_by_parts",
                          true,
                          "Whether to integrate the pressure term by parts.");
    params.addParam<bool>("supg",
                          false,
                          "Whether to perform SUPG stabilization of the momentum residuals");
    params.addParam<FunctionName>("forcing_func", 0, "The mms forcing function."););

template <ComputeStage compute_stage>
INSADMomentumBase<compute_stage>::INSADMomentumBase(const InputParameters & parameters)
  : INSADBase<compute_stage>(parameters),
    _component(adGetParam<unsigned>("component")),
    _integrate_p_by_parts(adGetParam<bool>("integrate_p_by_parts")),
    _supg(adGetParam<bool>("supg")),
    _ffn(this->getFunction("forcing_func"))
{
  if (_supg && !_convective_term)
    mooseError("It doesn't make sense to conduct SUPG stabilization without a convective term.");
}

template <ComputeStage compute_stage>
ADResidual
INSADMomentumBase<compute_stage>::computeQpResidual()
{
  ADResidual r = 0;

  // viscous term
  r += computeQpResidualViscousPart();

  // pressure term
  if (_integrate_p_by_parts)
    r += _grad_test[_i][_qp](_component) * this->weakPressureTerm();
  else
    r += _test[_i][_qp] * this->strongPressureTerm()(_component);

  // body force term
  r += _test[_i][_qp] * (this->gravityTerm()(_component) - _ffn.value(_t, _q_point[_qp]));

  // convective term
  if (_convective_term)
    r += _test[_i][_qp] * this->convectiveTerm()(_component);

  if (_supg)
    r += computeQpPGResidual();

  return r;
}

template <ComputeStage compute_stage>
ADResidual
INSADMomentumBase<compute_stage>::computeQpPGResidual()
{
  INSVectorValue<compute_stage> U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  const auto & convective_term =
      _convective_term ? this->convectiveTerm() : INSVectorValue<compute_stage>(0, 0, 0);
  const auto & viscous_term =
      _laplace ? this->strongViscousTermLaplace() : this->strongViscousTermTraction();
  const auto & transient_term =
      _transient_term ? this->timeDerivativeTerm() : INSVectorValue<compute_stage>(0, 0, 0);

  return this->tau() * U * _grad_test[_i][_qp] *
         ((convective_term + viscous_term + transient_term + this->strongPressureTerm() +
           this->gravityTerm())(_component)-_ffn.value(_t, _q_point[_qp]));
}

template class INSADMomentumBase<RESIDUAL>;
template class INSADMomentumBase<JACOBIAN>;
