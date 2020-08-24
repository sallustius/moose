#!/usr/bin/env python3

import mms
import sympy

# f_u, e_u = mms.evaluate('-div(mu * grad(u))', '1.1*sin(1.1*x)', variable='u', scalars=['mu'])
# f_u, e_u = mms.evaluate('div(vel*rho*u)', '1.1*sin(1.1*x)', variable='u', vel='1.1*sin(1.1*x) * e_i', scalars=['rho'])
f_u, e_u = mms.evaluate('div(vel*rho*u) - div(mu * grad(u))', '1.1*sin(1.1*x)', variable='u', vel='1.1*sin(1.1*x) * e_i', scalars=['mu', 'rho'])

rho = sympy.Symbol('rho')

e_rhou = e_u * rho

mms.print_hit(e_u, 'exact_u')
mms.print_hit(e_rhou, 'exact_rhou', rho='${rho}')
mms.print_hit(f_u, 'forcing_u', mu='${mu}', rho='${rho}')
