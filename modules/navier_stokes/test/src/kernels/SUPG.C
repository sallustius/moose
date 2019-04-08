//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SUPG.h"
#include "Function.h"

registerADMooseObject("NavierStokesTestApp", SUPG);

defineADValidParams(SUPG,
                    ADKernel,
                    params.addParam<FunctionName>("forcing_func",
                                                  0,
                                                  "The forcing function, typically used for MMS.");
                    params.addRequiredParam<RealVectorValue>("velocity", "The velocity");
                    params.addParam<Real>("diff", 0, "The diffusivity");
                    params.addParam<bool>("include_transient_term",
                                          false,
                                          "Whether to include the transient term in this kernel");
                    MooseEnum tau_type("none squares simple", "squares");
                    params.addRequiredParam<MooseEnum>(
                        "tau_type", tau_type, "The type of stabilization parameter to use."););

template <ComputeStage compute_stage>
SUPG<compute_stage>::SUPG(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _ffn(getFunction("forcing_func")),
    _velocity(adGetParam<RealVectorValue>("velocity")),
    _second_u(_var.template adSecondSln<compute_stage>()),
    _u_dot(_var.template adUDot<compute_stage>()),
    _diff(adGetParam<Real>("diff")),
    _include_transient_term(adGetParam<bool>("include_transient_term")),
    _tau_type(adGetParam<MooseEnum>("tau_type"))
{
}

template <ComputeStage compute_stage>
ADReal
SUPG<compute_stage>::computeQpResidual()
{
  auto h = _current_elem->hmax();
  ADReal tau;
  if (_tau_type == "squares")
    tau = 1. / std::sqrt(4. / (_dt * _dt) + 4. * _velocity * _velocity / (h * h) +
                         9. * (4. * _diff / (h * h)) * (4. * _diff / (h * h)));
  else if (_tau_type == "simple")
    tau = h / (2. * _velocity.norm());
  else if (_tau_type == "none")
    tau = 0;

  ADReal transient_term = _include_transient_term ? _u_dot[_qp] : 0;
  return tau * _velocity * _grad_test[_i][_qp] *
         (transient_term + _velocity * _grad_u[_qp] - _diff * _second_u[_qp].tr() -
          _ffn.value(_t, _q_point[_qp]));
}
