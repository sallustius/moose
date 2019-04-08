//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeDerivativeSUPG.h"
#include "Function.h"

registerADMooseObject("NavierStokesTestApp", TimeDerivativeSUPG);

defineADValidParams(TimeDerivativeSUPG,
                    ADKernel,
                    params.addRequiredParam<RealVectorValue>("velocity", "The velocity");
                    params.addParam<Real>("diff", 0, "The diffusivity");
                    params.set<MultiMooseEnum>("vector_tags") = "time";
                    params.set<MultiMooseEnum>("matrix_tags") = "system time";
                    MooseEnum tau_type("squares simple none", "squares");
                    params.addRequiredParam<MooseEnum>(
                        "tau_type", tau_type, "The type of stabilization parameter to use."););

template <ComputeStage compute_stage>
TimeDerivativeSUPG<compute_stage>::TimeDerivativeSUPG(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _velocity(adGetParam<RealVectorValue>("velocity")),
    _u_dot(_var.template adUDot<compute_stage>()),
    _diff(adGetParam<Real>("diff")),
    _tau_type(adGetParam<MooseEnum>("tau_type"))
{
}

template <ComputeStage compute_stage>
ADReal
TimeDerivativeSUPG<compute_stage>::computeQpResidual()
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

  return tau * _velocity * _grad_test[_i][_qp] * _u_dot[_qp];
}
