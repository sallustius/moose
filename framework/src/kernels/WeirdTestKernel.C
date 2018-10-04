//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeirdTestKernel.h"

registerMooseObject("MooseApp", WeirdTestKernel);

template <>
InputParameters
validParams<WeirdTestKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  params.addCoupledVar("disp_x", "The displacement");
  params.addCoupledVar("disp_y", "The displacement");
  params.addCoupledVar("disp_z", "The displacement");
  return params;
}

WeirdTestKernel::WeirdTestKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _disp_x_id(coupled("disp_x")),
    _disp_y_id(coupled("disp_y")),
    _disp_z_id(coupled("disp_z")),
    _dim(_mesh.dimension()),
    _dphidx_derivatives(_assembly.dphidxDerivatives()),
    _ad_grad_u(_var.adGradSln<JACOBIAN>())
{
}

Real
WeirdTestKernel::computeQpResidual()
{
  return 0;
}

void
WeirdTestKernel::computeResidual()
{
  prepareVectorTag(_assembly, _var.number());
  for (_i = 0; _i < _test.size(); _i++)
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      _local_re(_i) += _JxW[_qp] * _grad_u[_qp] * _grad_test[_i][_qp] * _coord[_qp];

  accumulateTaggedLocalResidual();
}

void
WeirdTestKernel::computeJacobian()
{
  prepareMatrixTag(_assembly, _var.number(), _var.number());

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        _local_ke(_i, _j) += _JxW[_qp] * _grad_phi[_j][_qp] * _grad_test[_i][_qp] * _coord[_qp];

  accumulateTaggedLocalMatrix();
}

void
WeirdTestKernel::computeOffDiagJacobian(MooseVariableFEBase & jvar)
{
  auto jvar_num = jvar.number();
  if (jvar_num == _var.number())
    computeJacobian();
  else
  {
    prepareMatrixTag(_assembly, _var.number(), jvar_num);

    if (_local_ke.m() != _test.size() || _local_ke.n() != jvar.phiSize())
      return;

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < jvar.phiSize(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
          if (jvar_num == _disp_x_id)
          {
            auto JxW_term =
                _JxW_derivatives[_qp][_dim * _j] * _grad_u[_qp] * _grad_test[_i][_qp] * _coord[_qp];
            auto grad_test_term = _JxW[_qp] * _grad_u[_qp] *
                                  RealVectorValue(_dphidx_derivatives[_i][_qp][_dim * _j], 0, 0) *
                                  _coord[_qp];
            auto grad_u_term = _JxW[_qp] *
                               RealVectorValue(_ad_grad_u[_qp].derivatives()[_dim * _j](0), 0, 0) *
                               _grad_test[_i][_qp] * _coord[_qp];
            _local_ke(_i, _j) += JxW_term + grad_test_term + grad_u_term;
          }
          else if (jvar_num == _disp_y_id)
            _local_ke(_i, _j) += 0;
          else if (jvar_num == _disp_y_id)
            _local_ke(_i, _j) += 0;
        }

    accumulateTaggedLocalMatrix();
  }
}
