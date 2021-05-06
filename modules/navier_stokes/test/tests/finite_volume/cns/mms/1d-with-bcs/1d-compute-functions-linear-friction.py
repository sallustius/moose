#!/usr/bin/env python3

import mms
import sympy

eps = 'cos(1.3*x)'

pressure = '1.01e5*cos(x)'
T = '273.15*cos(1.1*x)'
ud = '1.1*cos(1.2*x)'
gamma = 1.4
R = 8.3145
molar_mass = 29.0e-3
R_specific = R / molar_mass
cp = gamma * R_specific / (gamma - 1.)
cv = cp / gamma
rho = 'pressure * ' + str(molar_mass) + ' / ( ' + str(R) + ' * T)'
rho_ud = 'rho * ud'
u = 'ud / eps'
vel = 'u * e_i'
mass_flux = 'rho_ud * e_i'
e = 'pressure / (0.4 * rho)'
et = 'e + 0.5 * vel.dot(vel)'
ht = 'et + pressure / rho'

f_pressure, e_pressure = mms.evaluate('div(mass_flux)', pressure, variable='pressure', pressure=pressure, ud=ud, T=T, rho=rho, rho_ud=rho_ud, mass_flux=mass_flux)
f_ud, e_ud = mms.evaluate('div(mass_flux * u) + eps*grad(pressure).dot(e_i)', ud, variable='ud', pressure=pressure, ud=ud, T=T, eps=eps, rho=rho, rho_ud=rho_ud, mass_flux=mass_flux, u=u)
f_T, e_T = mms.evaluate('div(mass_flux * ht)', T, variable='T', pressure=pressure, ud=ud, T=T, eps=eps, rho=rho, rho_ud=rho_ud, mass_flux=mass_flux, e=e, u=u, vel=vel, et=et, ht=ht)

mms.print_hit(e_pressure, 'exact_p')
mms.print_hit(f_pressure, 'forcing_p')

mms.print_hit(e_ud, 'exact_ud')
mms.print_hit(f_ud, 'forcing_ud')

mms.print_hit(e_T, 'exact_T')
mms.print_hit(f_T, 'forcing_T')
