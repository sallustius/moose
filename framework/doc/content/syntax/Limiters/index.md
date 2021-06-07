# Limiters

Limiters, generally speaking, limit the slope when doing high-order (e.g. accuracy order greater than
1, e.g. non-constant polynomial) interpolations from finite volume cell
centroids to faces. This limiting is done to avoid creating oscillations in the
solution field in regions of steep gradients or discontinuities. Slope limiters,
or flux limiters, are generally employed to make the solution Total Variation
Diminishing (TVD). Borrowing notation from
[here](https://en.wikipedia.org/wiki/Total_variation_diminishing), the Total
Variation when space and time have been discretized can be defined as

\begin{equation}
TV(u^n) = TV(u(\centerdot,t^n)) = \sum_j \vert u_{j+1}^n - u_j^n \vert
\end{equation}

where $u$ is the discretized approximate solution, $n$ denotes the time index,
and $u_j^n = u(x_j,t^n)$. A numerical method is TVD if

\begin{equation}
TV(u^{n+1}) \leq TV(u^n)
\end{equation}


The following limiters are available in MOOSE. We have noted the convergence
orders of each (when considering that the solution is smooth) and whether they are TVD

Limiter class name | Convergence Order | TVD
------------------ | ----------------- | ---
`VanLeer`          | 2                 | Yes
`Upwind`           | 1                 | Yes
`CentralDifference` | 2                | No
`MinMod`           | 2                 | Yes
`SOU`              | 2                 | No
`QUICK`            | 2                 | No
