#!/usr/bin/env python3

import mms
import sympy

vel = '1.1*sin(1.1*x) * e_i'
p = '1.1*cos(1.1*x)'

f_u, e_u = mms.evaluate('div(vel*rho*u) - div(mu * grad(u)) + grad(p).dot(e_i)', '1.1*sin(1.1*x)', variable='u', vel=vel, p=p, scalars=['mu', 'rho'])

f_p, e_p = mms.evaluate('div(vel)', p, variable='p', vel=vel)

rho = sympy.Symbol('rho')

e_rhou = e_u * rho

mms.print_hit(e_u, 'exact_u')
mms.print_hit(e_rhou, 'exact_rhou', rho='${rho}')
mms.print_hit(f_u, 'forcing_u', mu='${mu}', rho='${rho}')
mms.print_hit(e_p, 'exact_p')
mms.print_hit(f_p, 'forcing_p')
