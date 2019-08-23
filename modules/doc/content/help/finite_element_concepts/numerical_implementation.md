# Numerical Integration

- The only remaining non-discretized parts of the weak form are the integrals.
- We split the domain integral into a sum of integrals over elements:
  $$$\int_{\Omega} f(\vec{x}) \;\text{d}\vec{x} = \sum_e \int_{\Omega_e} f(\vec{x}) \;\text{d}\vec{x}$$$
- Through a change of variables, the element integrals are mapped to integrals over the "reference" elements $$$\hat{\Omega}_e$$$.
  $$$\sum_e \int_{\Omega_e} f(\vec{x}) \;\text{d}\vec{x} =
        \sum_e \int_{\hat{\Omega}_e} f(\vec{\xi}) \left|\mathcal{J}_e\right| \;\text{d}\vec{\xi}$$$
- $$$\mathcal{J}_e$$$ is the Jacobian of the map from the physical element to the reference element.
- To approximate the reference element integrals numerically, we use quadrature (typically "Gaussian Quadrature"):
  $$$\sum_e \int_{\hat{\Omega}_e} f(\vec{\xi}) \left|\mathcal{J}_e\right| \;\text{d}\vec{\xi} \approx
        \sum_e \sum_{q} w_{q} f( \vec{x}_{q}) \left|\mathcal{J}_e(\vec{x}_{q})\right|$$$
- $$$\vec{x}_{q}$$$ is the spatial location of the $$$q$$$th quadrature point and $$$w_{q}$$$ is its associated weight.
- MOOSE handles multiplication by the Jacobian and the weight automatically, thus your `Kernel` is only responsible for computing the $$$f(\vec{x}_{q})$$$ part of the integrand.
- Under certain common situations, the quadrature approximation is exact!
    - For example, in 1 dimension, Gaussian Quadrature can exactly integrate polynomials of order $$$2n-1$$$ with $$$n$$$ quadrature points.
- Note that sampling $$$u_h$$$ at the quadrature points yields:
   $$$\begin{aligned}
    u(\vec{x}_{q}) &\approx u_h(\vec{x}_{q}) = \sum u_j \phi_j(\vec{x}_{q}) \\
    \nabla u (\vec{x}_{q}) &\approx \nabla u_h(\vec{x}_{q}) = \sum u_j \nabla \phi_j(\vec{x}_{q})\end{aligned}$$$
- And our weak form becomes:
  $$$\begin{aligned}
  R_i(u_h) &= \sum_e \sum_{q} w_{q} \left|\mathcal{J}_e\right|\left[ \nabla\psi_i\cdot k \nabla u_h + \psi_i \left(\vec\beta\cdot \nabla u_h \right) - \psi_i f \right](\vec{x}_{q}) \\
  &- \sum_f \sum_{q_{face}} w_{q_{face}} \left|\mathcal{J}_f\right| \left[\psi_i k \nabla u_h \cdot \vec{n} \right](\vec x_{q_{face}})
  \end{aligned}$$$
- The second sum is over boundary faces, $$$f$$$.
- MOOSE `Kernels` must provide each of the terms in square brackets (evaluated at $$$\vec{x}_{q}$$$ or $$$\vec x_{q_{face}}$$$ as necessary).

# Newton's Method

[NonlinearSystem.md#newtons_method]

# Newton for a Simple Equation

- Consider the convection-diffusion equation with nonlinear $$$k$$$, $$$\vec{\beta}$$$, and $$$f$$$:
    $$$\begin{aligned}- \nabla\cdot k\nabla u + \vec{\beta} \cdot \nabla u = f\end{aligned}$$$

- The $$$i^{th}$$$ component of the residual vector is:
    $$$\begin{aligned}
    R_i(u_h) = \left(\nabla\psi_i, k\nabla u_h \right) - \langle\psi_i, k\nabla u_h\cdot \hat{n} \rangle +
    \left(\psi_i, \vec{\beta} \cdot \nabla u_h\right) - \left(\psi_i, f\right)\end{aligned}$$$


- Using the previously-defined rules for $$$\frac{\partial u_h}{\partial u_j}$$$ and $$$\frac{\partial \left(\nabla u_h\right)}{\partial u_j}$$$, the $$$(i,j)$$$ entry of the Jacobian is then:

$$$\begin{aligned} J_{ij}(u_h) &= \left(\nabla\psi_i, \frac{\partial k}{\partial u_j}\nabla u_h \right) + \left(\nabla\psi_i, k \nabla \phi_j \right) - \left \langle\psi_i, \frac{\partial k}{\partial u_j}\nabla u_h\cdot \hat{n} \right\rangle \\&- \left \langle\psi_i, k\nabla \phi_j\cdot \hat{n} \right\rangle + \left(\psi_i, \frac{\partial \vec{\beta}}{\partial u_j} \cdot\nabla u_h\right) + \left(\psi_i, \vec{\beta} \cdot \nabla \phi_j\right) - \left(\psi_i, \frac{\partial f}{\partial u_j}\right)\end{aligned}$$$

- Note that even for this "simple" equation, the Jacobian entries are nontrivial: they depend on the partial derivatives of $$$k$$$, $$$\vec{\beta}$$$, and $$$f$$$, which may be difficult or time-consuming to compute analytically.

- In a multiphysics setting with many coupled equations and complicated material properties, the Jacobian might be extremely difficult to determine.

# Chain Rule

- On the previous slide, the term $$$\frac{\partial f}{\partial u_j}$$$ was used, where $$$f$$$ was a nonlinear forcing function.

- The chain rule allows us to write this term as

  $$$\begin{aligned}
    \frac{\partial f}{\partial u_j} &= \frac{\partial f}{\partial u_h} \frac{\partial u_h}{\partial u_j}
    \\
    &=\frac{\partial f}{\partial u_h} \phi_j\end{aligned}$$$

- If a functional form of $$$f$$$ is known, e.g. $$$f(u) = \sin(u)$$$, this
  formula implies that its Jacobian contribution is given by

   $$$\frac{\partial f}{\partial u_j} = \cos(u_h) \phi_j$$$

# Jacobian Free Newton Krylov

[NonlinearSystem.md#JFNK]

# Wrap Up

- The Finite Element Method is a way of numerically approximating the solution of PDEs.
- Just like polynomial fitting, FEM finds coefficients for basis functions.
- The "solution" is the combination of the coefficients and the basis functions, and the solution can be sampled anywhere in the domain.
- We compute integrals numerically using quadrature.
- Newton's Method provides a mechanism for solving a system of nonlinear equations.
- The Jacobian Free Newton Krylov (JFNK) method allows us to avoid explicitly forming the Jacobian matrix while still computing its "action".
