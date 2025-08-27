# Numerics

## Temporal discretization

### Pseudo-incompressible mode

The pseudo-incompressible equations are integrated with the semi-implicit time scheme

$$\begin{align*}
    1. && \rho'^n & = \rho^n - \overline{\rho},\\
    2. && \left(\rho^\#, \rho'^\#, \widehat{\boldsymbol{u}}^\#\right) & = \mathrm{L}_{\Delta t / 2} \left(\rho^n, \rho'^n, \widehat{\boldsymbol{u}}^n, \widehat{\boldsymbol{u}}^n, \alpha_\mathrm{R}^{n + 1}\right)\\
    3. && \left(\rho'^{n + 1 / 2}, \widehat{\boldsymbol{u}}^{n + 1 / 2}, \pi'^{n + 1 / 2}\right) & = \mathrm{RI}_{\Delta t / 2} \left(\rho^\#, \rho'^\#, \widehat{\boldsymbol{u}}^\#, \pi'^n, \beta_\mathrm{R}^{uv, n + 1}, \beta_\mathrm{R}^{\widehat{w}, n + 1}\right)\\
    4. && \left(\rho'^*, \widehat{\boldsymbol{u}}^*\right) & = \mathrm{RE}_{\Delta t / 2} \left(\rho^n, \rho'^n, \widehat{\boldsymbol{u}}^n, \pi'^{n + 1 / 2}\right)\\
    5. && \left(\rho^{**}, \rho'^{**}, \widehat{\boldsymbol{u}}^{**}\right) & = \mathrm{L}_{\Delta t} \left(\rho^n, \rho'^*, \widehat{\boldsymbol{u}}^*, \widehat{\boldsymbol{u}}^{n + 1 / 2}, \alpha_\mathrm{R}^{n + 1}\right)\\
    6. && \left(\rho'^{n + 1}, \widehat{\boldsymbol{u}}^{n + 1}, \pi'^{n + 1}\right) & = \mathrm{RI}_{\Delta t / 2} \left(\rho^{**}, \rho'^{**}, \widehat{\boldsymbol{u}}^{**}, \pi'^{n + 1 / 2}, 2 \beta_\mathrm{R}^{uv, n + 1}, 2 \beta_\mathrm{R}^{\widehat{w}, n + 1}\right),
\end{align*}$$

where the operators $\mathrm{L}$, $\mathrm{RI}$ and $\mathrm{RE}$ perform an explicit integration of the left-hand sides, an implicit integration of the right-hand sides and an explicit integration of the right-hand sides, each over the time step indicated in its subscript, respectively. The superscripts represent various time levels between those before ($n$) and after ($n + 1$) the current time step $\Delta t = t^{n + 1} - t^n$. In $\mathrm{L}$, the fourth argument is the velocity by which the prognostic variables are transported ([Benacchio & Klein, 2019](https://doi.org/10.1175/mwr-d-19-0073.1); [Schmid et al., 2021](https://doi.org/10.1175/MWR-D-21-0126.1)). A complete description of the exact implementation of these steps follows here.

 1. The density fluctuations are synchronized with the full density, i.e.

    $$\rho'^n = \rho^n - \overline{\rho}.$$

 1. The left-hand sides are integrated over $\Delta t / 2$ with a low-storage RK3 scheme ([Williamson, 1980](https://doi.org/10.1016/0021-9991(80)90033-9)). Fractional implicit Euler steps are used to integrate the Rayleigh-damping terms of the LHS sponge. At every RK3 stage $m$, the following updates are performed.

     1. Density update:

        $$\begin{align*}
            q^{\rho, m + 1} & = - \frac{\Delta t}{2 J} \left(\frac{\partial J \rho^m u^n}{\partial \widehat{x}} + \frac{\partial J \rho^m v^n}{\partial \widehat{y}} + \frac{\partial J \rho^m \widehat{w}^n}{\partial \widehat{z}}\right) + \left(\alpha_\mathrm{RK} q^\rho\right)^m,\\
            \rho^{m + 1} & = \rho^m + \beta_\mathrm{RK}^m q^{\rho, m + 1},\\
            \rho^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \left(\rho^{m + 1} + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \overline{\rho}\right)
        \end{align*}$$

     1. Density-fluctuations update:

        $$\begin{align*}
            q^{\rho', m + 1} & = - \frac{\Delta t}{2 J} \left(\frac{\partial J \rho'^m u^n}{\partial \widehat{x}} + \frac{\partial J \rho'^m v^n}{\partial \widehat{y}} + \frac{\partial J \rho'^m \widehat{w}^n}{\partial \widehat{z}}\right) + \left(\alpha_\mathrm{RK} q^{\rho'}\right)^m,\\
            \rho'^{m + 1} & = \rho'^m + \beta_\mathrm{RK}^m q^{\rho', m + 1},\\
            \rho'^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \rho'^{m + 1}
        \end{align*}$$

     1. Wind update:

        $$\begin{align*}
            q^{\rho u, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{u, m, n} + \mathcal{V}^{u, m} + f \left(\rho v\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho u}\right)^m,\\
            u^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho u\right)^m + \beta_\mathrm{RK}^m q^{\rho u, m + 1}\right],\\
            q^{\rho v, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{v, m, n} + \mathcal{V}^{v, m} - f \left(\rho u\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho v}\right)^m,\\
            v^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho v\right)^m + \beta_\mathrm{RK}^m q^{\rho v, m + 1}\right],\\
            q^{\rho \widehat{w}, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{\widehat{w}, m, n} + \mathcal{V}^{\widehat{w}, m} + G^{13} f \left(\rho v\right)^m - G^{23} f \left(\rho u\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho \widehat{w}}\right)^m\\
            \widehat{w}^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho \widehat{w}\right)^m + \beta_\mathrm{RK}^m q^{\rho \widehat{w}, m + 1}\right],\\
            \widehat{\boldsymbol{u}}^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \left(\widehat{\boldsymbol{u}}^{m + 1} + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \widehat{\boldsymbol{u}}_\mathrm{R}^{m + 1}\right)
        \end{align*}$$

    Therein, $\alpha_\mathrm{RK}^m$ and $\beta_\mathrm{RK}^m$ are the RK3 coefficients, with $\alpha_\mathrm{RK}^m = 0$ at the first stage, and $f_\mathrm{RK}^m$ are the ratios between the sizes of the full time step and the RK3 substeps. The advective momentum-flux divergences are given by

    $$\begin{align*}
        \mathcal{A}^{u, m, n} & = \frac{1}{J} \left[\frac{\partial J \left(\rho u\right)^m u^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho u\right)^m v^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho u\right)^m \widehat{w}^n}{\partial \widehat{z}}\right],\\
        \mathcal{A}^{v, m, n} & = \frac{1}{J} \left[\frac{\partial J \left(\rho v\right)^m u^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho v\right)^m v^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho v\right)^m \widehat{w}^n}{\partial \widehat{z}}\right],\\
        \mathcal{A}^{\widehat{w}, m, n} & = G^{13} \mathcal{A}^{u, m, n} + G^{23} \mathcal{A}^{v, m, n} + \frac{1}{J^2} \left[\frac{\partial J \left(\rho w\right)^m u^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho w\right)^m v^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho w\right)^m \widehat{w}^n}{\partial \widehat{z}}\right].
    \end{align*}$$

 1. The right-hand sides are integrated over $\Delta t / 2$ with an implicit Euler step. The divergence constraint is then enforced by solving a Poisson equation for Exner-pressure differences $\Delta \pi'$, which are used to correct the wind, the the density fluctuations and the Exner-pressure itself. The details of these substeps are as follows.

     1. Predictor step:

        $$\begin{align*}
            u^{n + 1 / 2} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \left[u^\# + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^\#} \mathcal{P}^{u, n} + \frac{F^{u, n + 1}}{\rho^\#}\right)\right],\\
            v^{n + 1 / 2} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \left[v^\# + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^\#} \mathcal{P}^{v, n} + \frac{F^{v, n + 1}}{\rho^\#}\right)\right],\\
            \widehat{w}^{n + 1 / 2} & = \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left[\widehat{w}^\# + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^\#} \mathcal{P}^{\widehat{w}, n} + \frac{b'^\#}{J} + \frac{F^{\widehat{w}, n + 1}}{\rho^\#}\right) \vphantom{\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}}\right.\\
            & \qquad \quad + \left.\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4} \left(G^{1 3} u^{n + 1 / 2} + G^{2 3} v^{n + 1 / 2}\right)\right],\\
            \rho'^{n + 1 / 2} & = - \frac{\rho^\#}{g} \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left\{- \frac{\overline{\rho}}{\rho^\#} N^2 \frac{\Delta t}{2} J \left[\widehat{w}^\# + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^\#} \mathcal{P}^{\widehat{w}, n} + \frac{F^{\widehat{w}, n + 1}}{\rho^\#}\right)\right]\right.\\
            & \qquad \quad + \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right) b'^\# + \frac{\overline{\rho}}{\rho^\#} N^2 \frac{\Delta t}{2} J \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right)\\
            & \qquad \quad \times \left.\left(G^{1 3} u^{n + 1 / 2} + G^{2 3} v^{n + 1 / 2}\right) \vphantom{\frac{\overline{\rho}}{\rho^\#} N^2 \frac{\Delta t}{2}}\right\}
        \end{align*}$$

     1. Poisson equation:

        $$\begin{align*}
            & \frac{1}{J} \left(\frac{\partial J P u^{n + 1 / 2}}{\partial \widehat{x}} + \frac{\partial J P v^{n + 1 / 2}}{\partial \widehat{y}} + \frac{\partial J P \widehat{w}^{n + 1 / 2}}{\partial \widehat{z}}\right)\\
            & \quad = \frac{1}{J} \left\{\frac{\partial J P C^{u, \#, n + 1}}{\partial \widehat{x}} + \frac{\partial J P C^{v, \#, n + 1}}{\partial \widehat{y}} \vphantom{\left(\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right)^{- 1}}\right.\\
            & \qquad \qquad + \frac{\partial}{\partial \widehat{z}} \left[\left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right)^{- 1}\right.\\
            & \qquad \qquad \qquad \quad \times \left(\frac{\Delta t}{2} c_p \frac{J P^2}{\rho^\#} \mathcal{D}^{\widehat{w}} + J P \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right.\\
            & \qquad \qquad \qquad \qquad \quad \times \left.\left.\left.\left(G^{13} \mathcal{C}^{u, \#, n + 1} + G^{23} \mathcal{C}^{v, \#, n + 1}\right) \vphantom{\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}}\right) \vphantom{\left(\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right)^{- 1}}\right]\right\}
        \end{align*}$$

     1. Corrector step:

        $$\begin{align*}
            u^{n + 1 / 2} & \rightarrow u^{n + 1 / 2} - \mathcal{C}^{u, \#, n + 1},\\
            v^{n + 1 / 2} & \rightarrow v^{n + 1 / 2} - \mathcal{C}^{v, \#, n + 1},\\
            \widehat{w}^{n + 1 / 2} & \rightarrow \widehat{w}^{n + 1 / 2} - \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left[\frac{\Delta t}{2} c_p \frac{P}{\rho^\#} \mathcal{D}^{\widehat{w}} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4} \left(G^{1 3} \mathcal{C}^{u, \#, n + 1} + G^{23} \mathcal{C}^{v, \#, n + 1}\right)\right],\\
            \rho'^{n + 1 / 2} & \rightarrow \rho'^{n + 1 / 2} + \frac{\rho^\#}{g} \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left[- \frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4} J c_p \frac{P}{\rho^\#} \mathcal{D}^{\widehat{w}} + \frac{\overline{\rho}}{\rho^\#} N^2 \frac{\Delta t}{2} J \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right)\right.\\
            & \qquad \quad \times \left.\left(G^{1 3} \mathcal{C}^{u, \#, n + 1} + G^{2 3} \mathcal{C}^{v, \#, n + 1}\right) \vphantom{\frac{\overline{\rho}}{\rho^\#} \frac{\left(N \Delta t\right)^2}{4}}\right],\\
            \pi'^{n + 1 / 2} & = \pi'^n + \Delta \pi'
        \end{align*}$$

    Therein, $b' = - g \rho' / \rho$, the pressure-difference gradient is given by

    $$\begin{align*}
        \mathcal{D}^u & = \left(\frac{\partial \Delta \pi'}{\partial \widehat{x}} + G^{13} \frac{\partial \Delta \pi'}{\partial \widehat{z}}\right),\\
        \mathcal{D}^v & = \left(\frac{\partial \Delta \pi'}{\partial \widehat{y}} + G^{23} \frac{\partial \Delta \pi'}{\partial \widehat{z}}\right),\\
        \mathcal{D}^{\widehat{w}} & = \left(G^{13} \frac{\partial \Delta \pi'}{\partial \widehat{x}} + G^{23} \frac{\partial \Delta \pi'}{\partial \widehat{y}} + G^{33} \frac{\partial \Delta \pi'}{\partial \widehat{z}}\right)
    \end{align*}$$

    and the horizontal-wind correction terms are

    $$\begin{align*}
        \mathcal{C}^{u, \#, n + 1} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \frac{\Delta t}{2} c_p \frac{P}{\rho^\#} \mathcal{D}^u,\\
        \mathcal{C}^{v, \#, n + 1} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \frac{\Delta t}{2} c_p \frac{P}{\rho^\#} \mathcal{D}^v.
    \end{align*}$$

 1. The right-hand sides (without the Rayleigh-damping terms) are integrated over $\Delta t / 2$ with an explicit Euler step. Therein, the Exner-pressure fluctuations $\pi'^{n + 1 / 2}$ are used to compute the pressure gradient. The exact updates are

    $$\begin{align*}
        u^* & = u^n + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^n} \mathcal{P}^{u, n + 1 / 2} + \frac{F^{u, n + 1}}{\rho^n}\right),\\
        v^* & = v^n + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^n} \mathcal{P}^{v, n + 1 / 2} + \frac{F^{v, n + 1}}{\rho^n}\right),\\
        \widehat{w}^* & = \widehat{w}^n + \frac{\Delta t}{2} \left(- c_p \frac{P}{\rho^n} \mathcal{P}^{\widehat{w}, n + 1 / 2} + \frac{b'^n}{J} + \frac{F^{\widehat{w}, n + 1}}{\rho^n}\right),\\
        \rho'^* & = - \frac{\rho^n}{g} \left[b'^n - \frac{\Delta t}{2} \overline{\rho} N^2 \left(\frac{w}{\rho}\right)^n\right].
    \end{align*}$$

 1. The left-hand sides are integrated over $\Delta t$ with the low-storage RK3 scheme. Fractional implicit Euler steps are once again used to integrate the Rayleigh-damping terms of the LHS sponge. This step is equivalent to the first one, except for the differences indicated in the compact description above.

 1. The right-hand sides are integrated over $\Delta t / 2$ with an implicit Euler step, followed by the Poisson equation being solved and a correction step being performed. The Rayleigh-damping terms are doubled, since they were left out in the explicit Euler step. This step is equivalent to the second one, except for the differences indicated in the compact description above.

### Boussinesq mode

In Boussinesq mode, the time scheme remains mostly unchanged. As has been mentioned in the description of the physics, the continuity equation is removed, as are the density fluctuations everywhere except in the auxiliary equation and the buoyancy term of the transformed-vertical-momentum equation. Furthermore, $\overline{\rho}$, $\overline{\theta}$, $P$ and $N^2$ are replaced with $\rho_0$, $\theta_0$, $P_0$ and $N_0^2$, respectively.

### Compressible mode

In compressible mode, the time scheme is changed in several ways, due to the mass-weighted potential temperature having a spatiotemporal dependence and being directly coupled to the Exner-pressure. It may be summarized by

$$\begin{align*}
    1. && \rho'^n & = \rho^n - \frac{P^n}{\overline{\theta}},\\
    2. && \left(\rho^\#, \rho'^\#, P^\#, \widehat{\boldsymbol{u}}^\#, \pi'^\#\right) & = \mathrm{L}_{\Delta t / 2} \left(\rho^n, \rho'^n, P^n, \widehat{\boldsymbol{u}}^n, \pi'^n, P^n, \widehat{\boldsymbol{u}}^n, \alpha_\mathrm{R}^{n + 1}\right)\\
    3. && \left(\rho'^{n + 1 / 2}, \widehat{\boldsymbol{u}}^{n + 1 / 2}, \pi'^{n + 1 / 2}\right) & = \mathrm{RI}_{\Delta t / 2} \left(\rho^\#, \rho'^\#, P^\#, \widehat{\boldsymbol{u}}^\#, \pi'^\#, \beta_\mathrm{R}^{uv, n + 1}, \beta_\mathrm{R}^{\widehat{w}, n + 1}\right)\\
    4. && \left(\rho'^*, \widehat{\boldsymbol{u}}^*, \pi'^*\right) & = \mathrm{RE}_{\Delta t / 2} \left(\rho^n, \rho'^n, P^n, \widehat{\boldsymbol{u}}^n, \pi'^n\right)\\
    5. && \left(\rho^{**}, \rho'^{**}, P^{**}, \widehat{\boldsymbol{u}}^{**}, \pi'^{**}\right) & = \mathrm{L}_{\Delta t} \left(\rho^n, \rho'^*, P^n, \widehat{\boldsymbol{u}}^*, \pi'^{n + 1 / 2}, P^\#, \widehat{\boldsymbol{u}}^{n + 1 / 2}, \alpha_\mathrm{R}^{n + 1}\right)\\
    6. && \left(\rho'^{n + 1}, \widehat{\boldsymbol{u}}^{n + 1}, \pi'^{n + 1}\right) & = \mathrm{RI}_{\Delta t / 2} \left(\rho^{**}, \rho'^{**}, P^{**}, \widehat{\boldsymbol{u}}^{**}, \pi'^{**}, 2 \beta_\mathrm{R}^{uv, n + 1}, 2 \beta_\mathrm{R}^{\widehat{w}, n + 1}\right)
\end{align*}$$

(see [Benacchio & Klein, 2019](https://doi.org/10.1175/mwr-d-19-0073.1); [Chew et al., 2022](https://doi.org/10.1175/MWR-D-21-0175.1)). Another detailed description follows here.

 1. The density fluctuations are synchronized with the full density, i.e.

    $$\rho'^n = \rho^n - \frac{P^n}{\overline{\theta}}.$$

 1. The left-hand sides are integrated over $\Delta t / 2$ with the low-storage RK3 scheme. Fractional implicit Euler steps are used to integrate the Rayleigh-damping terms of the LHS sponge. The fact that this is performed such that $P \widehat{\boldsymbol{u}}$ is the transporting wind now becomes relevant not only for the spatial discretization but also for the temporal one, so that the updates at every RK3 stage are as follows.

     1. Density update:

        $$\begin{align*}
            q^{\rho, m + 1} & = - \frac{\Delta t}{2 J} \left[\frac{\partial J \left(\rho / P\right)^m \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho / P\right)^m \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho / P\right)^m \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right] + \left(\alpha_\mathrm{RK} q^\rho\right)^m,\\
            \rho^{m + 1} & = \rho^m + \beta_\mathrm{RK}^m q^{\rho, m + 1},\\
            \rho^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \left(\rho^{m + 1} + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \overline{\rho}\right)
        \end{align*}$$

     1. Density-fluctuations update:

        $$\begin{align*}
            q^{\rho', m + 1} & = - \frac{\Delta t}{2} \left\{\frac{1}{J} \left[\frac{\partial J \left(\rho' / P\right)^m \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho' / P\right)^m \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho' / P\right)^m \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right]\right.\\
            & \qquad \qquad + \left.\frac{F^{P, n + 1}}{\overline{\theta}}\right\} + \left(\alpha_\mathrm{RK} q^{\rho'}\right)^m,\\
            \rho'^{m + 1} & = \rho'^m + \beta_\mathrm{RK}^m q^{\rho', m + 1},\\
            \rho'^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \left\{\rho'^{m + 1} + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \overline{\rho} \left[1 - \frac{\left(P / \rho\right)^{m + 1}}{\overline{\theta}}\right]\right\}
        \end{align*}$$

     1. Mass-weighted potential-temperature update:

        $$\begin{align*}
            q^{P, m + 1} & = - \frac{\Delta t}{2} \left\{\frac{1}{J} \left[\frac{\partial J \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right] - F^{P, n + 1}\right\}\\
            & \quad + \left(\alpha_\mathrm{RK} q^P\right)^m,\\
            P^{m + 1} & = P^m + \beta_\mathrm{RK}^m q^{P, m + 1},\\
            P^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} P^{m + 1} \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \frac{\overline{\rho}}{\rho^{m + 1}}\right)
        \end{align*}$$

     1. Wind update:

        $$\begin{align*}
            q^{\rho u, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{u, m, n} + \mathcal{V}^{u, m} + f \left(\rho v\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho u}\right)^m,\\
            u^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho u\right)^m + \beta_\mathrm{RK}^m q^{\rho u, m + 1}\right],\\
            q^{\rho v, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{v, m, n} + \mathcal{V}^{v, m} - f \left(\rho u\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho v}\right)^m,\\
            v^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho v\right)^m + \beta_\mathrm{RK}^m q^{\rho v, m + 1}\right],\\
            q^{\rho \widehat{w}, m + 1} & = \frac{\Delta t}{2} \left[- \mathcal{A}^{\widehat{w}, m, n} + \mathcal{V}^{\widehat{w}, m} + G^{13} f \left(\rho v\right)^m - G^{23} f \left(\rho u\right)^m\right] + \left(\alpha_\mathrm{RK} q^{\rho \widehat{w}}\right)^m,\\
            \widehat{w}^{m + 1} & = \left(\rho^{- 1}\right)^{m + 1} \left[\left(\rho \widehat{w}\right)^m + \beta_\mathrm{RK}^m q^{\rho \widehat{w}, m + 1}\right],\\
            \widehat{\boldsymbol{u}}^{m + 1} & \rightarrow \left(1 + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2}\right)^{- 1} \left(\widehat{\boldsymbol{u}}^{m + 1} + \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \frac{\Delta t}{2} \widehat{\boldsymbol{u}}_\mathrm{R}^{m + 1}\right)
        \end{align*}$$

    The advective momentum-flux divergences are given by

    $$\begin{align*}
        \mathcal{A}^{u, m, n} & = \frac{1}{J} \left[\frac{\partial J \left(\rho u / P\right)^m \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho u / P\right)^m \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho u / P\right)^m \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right],\\
        \mathcal{A}^{v, m, n} & = \frac{1}{J} \left[\frac{\partial J \left(\rho v / P\right)^m \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho v / P\right)^m \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho v / P\right)^m \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right],\\
        \mathcal{A}^{\widehat{w}, m, n} & = G^{13} \mathcal{A}^{u, m, n} + G^{23} \mathcal{A}^{v, m, n}\\
        & \quad + \frac{1}{J^2} \left[\frac{\partial J \left(\rho w / P\right)^m \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(\rho w / P\right)^m \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(\rho w / P\right)^m \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right].
    \end{align*}$$

    After the last RK3 step, the Exner-pressure fluctuations are updated according to the Rayleigh damping applied to the mass-weighted potential temperature, i.e.

    $$\begin{align*}
        \pi'^\# = \pi'^n - \alpha_\mathrm{R} \frac{\Delta t}{2} \left(P \frac{\partial \pi'}{\partial P}\right)^\# \left(1 - \frac{\overline{\rho}}{\rho^\#}\right),
    \end{align*}$$

    where $\left(\partial P / \partial \pi'\right)^\# = \left(\gamma - 1\right)^{- 1} \left(R / p_\mathrm{ref}\right)^{1 - \gamma} \left(P^{2 - \gamma}\right)^\#$, with $\gamma$ being the ratio between the specific heat capacities and constant pressure and volume, $R$ being the specific gas constant and $p_\mathrm{ref}$ being the reference (ground) pressure.

 1. The right-hand sides are integrated over $\Delta t / 2$ with an implicit Euler step. The divergence constraint is then enforced by solving a Poisson equation for Exner-pressure differences $\Delta \pi'$, which are used to correct the wind, the density fluctuations and the Exner-pressure itself. In principle, this is the same as in pseudo-incompressible mode, however, instead of $\widehat{\boldsymbol{u}}^{n + 1 / 2}$, the predictor step calculates $\widehat{\boldsymbol{U}}^{n + 1 / 2}$, which is obtained by multiplying the equation with $P^\#$. The details are as follows.

     1. Predictor step:

        $$\begin{align*}
            U^{n + 1 / 2} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \left\{\left(P u\right)^\# + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^\# \mathcal{P}^{u, n} + \left(\frac{P}{\rho}\right)^\# F^{u, n + 1}\right]\right\},\\
            V^{n + 1 / 2} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \left\{\left(P v\right)^\# + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^\# \mathcal{P}^{v, n} + \left(\frac{P}{\rho}\right)^\# F^{v, n + 1}\right]\right\},\\
            \widehat{W}^{n + 1 / 2} & = \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left\{\left(P \widehat{w}\right)^\# + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^\# \mathcal{P}^{\widehat{w}, n} + \frac{\left(P b'\right)^\#}{J} + \left(\frac{P}{\rho}\right)^\# F^{\widehat{w}, n + 1}\right]\right.\\
            & \qquad \quad + \left.\frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4} \left[G^{1 3} U^{n + 1 / 2} + G^{2 3} V^{n + 1 / 2}\right]\right\},\\
            \rho'^{n + 1 / 2} & = - \frac{\rho^\#}{g} \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left\{- \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} N^2 \frac{\Delta t}{2} J \left[\widehat{w}^\# + \frac{\Delta t}{2} \left(- c_p \left(\frac{P}{\rho}\right)^\# \mathcal{P}^{\widehat{w}, n} + \frac{F^{\widehat{w}, n + 1}}{\rho^\#}\right)\right]\right.\\
            & \qquad \quad  + \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right) b'^\# + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} N^2 \frac{\Delta t}{2} J \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right)\\
            & \qquad \quad \times \left.\left(G^{1 3} \frac{U^{n + 1 / 2}}{P^\#} + G^{2 3} \frac{V^{n + 1 / 2}}{P^\#}\right) \vphantom{\left[\left(\left(\frac{P}{\rho}\right)^\#\right)\right]}\right\}
        \end{align*}$$

     1. Poisson equation:

        $$\begin{align*}
            & \frac{1}{J} \left[\frac{\partial J U^{n + 1 / 2}}{\partial \widehat{x}} + \frac{\partial J V^{n + 1 / 2}}{\partial \widehat{y}} + \frac{\partial J \widehat{W}^{n + 1 / 2}}{\partial \widehat{z}}\right] - F^{P, n + 1}\\
            & \quad = - \left(\frac{\partial P}{\partial \pi'}\right)^\# \frac{2 \Delta \pi'}{\Delta t}\\
            & \qquad + \frac{1}{J} \left\{\frac{\partial J P^\# C^{u, \#, n + 1}}{\partial \widehat{x}} + \frac{\partial J P^\# C^{v, \#, n + 1}}{\partial \widehat{y}} \vphantom{\left[\left(\frac{\left(P /  \rho\right)^\#}{\overline{\theta}}\right)^{- 1}\right]}\right.\\
            & \qquad \qquad \quad + \frac{\partial}{\partial \widehat{z}} \left[\left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right)^{- 1}\right.\\
            & \qquad \qquad \qquad \qquad \times \left(\frac{\Delta t}{2} J c_p \left(\frac{P^2}{\rho}\right)^\# \mathcal{D}^{\widehat{w}} + J \frac{\left(P^2 /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right.\\
            & \qquad \qquad \qquad \qquad \qquad \times \left.\left.\left.\left(G^{13} \mathcal{C}^{u, \#, n + 1} + G^{23} \mathcal{C}^{v, \#, n + 1}\right)\vphantom{\left(\frac{P^2}{\rho}\right)^\#}\right) \vphantom{\left(\frac{\left(P /  \rho\right)^\#}{\overline{\theta}}\right)^{- 1}}\right]\right\}
        \end{align*}$$

     1. Corrector step:

        $$\begin{align*}
            u^{n + 1 / 2} & = \left(P^{- 1}\right)^\# \left(U^{n + 1 / 2} - P^\# \mathcal{C}^{u, \#, n + 1}\right),\\
            v^{n + 1 / 2} & = \left(P^{- 1}\right)^\# \left(V^{n + 1 / 2} - P^\# \mathcal{C}^{v, \#, n + 1}\right),\\
            \widehat{w}^{n + 1 / 2} & = \left(P^{- 1}\right)^\# \left\{\widehat{W}^{n + 1 / 2} - \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\right.\\
            & \qquad \qquad \qquad \times \left[\frac{\Delta t}{2} c_p \left(\frac{P^2}{\rho}\right)^\# \mathcal{D}^{\widehat{w}} + \frac{\left(P^2 /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right.\\
            & \qquad \qquad \qquad \qquad \times \left.\left.\left(G^{1 3} \mathcal{C}^{u, \#, n + 1} + G^{23} \mathcal{C}^{v, \#, n + 1}\right) \vphantom{\left(\frac{P^2}{\rho}\right)^\#}\right] \vphantom{\left[\frac{\left(P /  \rho\right)^\#}{\overline{\theta}}\right]^{- 1}}\right\},\\
            \rho'^{n + 1 / 2} & \rightarrow \rho'^{n + 1 / 2} + \frac{\rho^\#}{g} \left[1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2} + \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4}\right]^{- 1}\\
            & \quad \times \left[- \frac{\left(P /  \rho\right)^\#}{\overline{\theta}} \frac{\left(N \Delta t\right)^2}{4} J c_p \left(\frac{P}{\rho}\right)^\# \mathcal{D}^{\widehat{w}}\right.\\
            & \qquad \quad + \left.\frac{\left(P /  \rho\right)^\#}{\overline{\theta}} N^2 \frac{\Delta t}{2} J \left(1 + \beta_\mathrm{R}^{\widehat{w}, n + 1} \frac{\Delta t}{2}\right)  \left(G^{1 3} \mathcal{C}^{u, \#, n + 1} + G^{2 3} \mathcal{C}^{v, \#, n + 1}\right)\right],\\
            \pi'^{n + 1 / 2} & = \pi'^n + \Delta \pi'
        \end{align*}$$

    The horizontal-wind correction terms are given by

    $$\begin{align*}
        \mathcal{C}^{u, \#, n + 1} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \frac{\Delta t}{2} c_p \left(\frac{P}{\rho}\right)^\# \mathcal{D}^u,\\
        \mathcal{C}^{v, \#, n + 1} & = \left(1 + \beta_\mathrm{R}^{uv, n + 1} \frac{\Delta t}{2}\right)^{- 1} \frac{\Delta t}{2} c_p \left(\frac{P}{\rho}\right)^\# \mathcal{D}^v.
    \end{align*}$$

 1. The right-hand sides (without the Rayleigh-damping terms) are integrated over $\Delta t / 2$ with an explicit Euler step. In contrast to the pseudo-incompressible mode, the Exner-pressure fluctuations $\pi'^n$ are used to compute the pressure gradient. They are themselves also updated with the left-hand side of the prognostic equation for $P$. Thus, one has

    $$\begin{align*}
        u^* & = \left(P^{- 1}\right)^n \left\{\left(P u\right)^n + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^n \mathcal{P}^{u, n} + \left(\frac{P}{\rho}\right)^n F^{u, n + 1}\right]\right\},\\
        v^* & = \left(P^{- 1}\right)^n \left\{\left(P v\right)^n + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^n \mathcal{P}^{v, n} + \left(\frac{P}{\rho}\right)^n F^{v, n + 1}\right]\right\},\\
        \widehat{w}^* & = \left(P^{- 1}\right)^n \left\{\left(P \widehat{w}\right)^n + \frac{\Delta t}{2} \left[- c_p \left(\frac{P^2}{\rho}\right)^n \mathcal{P}^{\widehat{w}, n} + \frac{\left(P b'\right)^n}{J} + \left(\frac{P}{\rho}\right)^n F^{\widehat{w}, n + 1}\right]\right\},\\
        \rho'^* & = - \frac{\rho^n}{g} \left[b'^n - \frac{\Delta t}{2} N^2 \frac{\left(P w / \rho\right)^n}{\overline{\theta}}\right],\\
        \pi'^* & = \pi'^n - \frac{\Delta t}{2} \left(\frac{\partial \pi'}{\partial P}\right)^n \left\{\frac{1}{J} \left[\frac{\partial J \left(P u\right)^n}{\partial \widehat{x}} + \frac{\partial J \left(P v\right)^n}{\partial \widehat{y}} + \frac{\partial J \left(P \widehat{w}\right)^n}{\partial \widehat{z}}\right] - F^{P, n + 1}\right\}.
    \end{align*}$$

 1. The left-hand sides are integrated over $\Delta t$ with the low-storage RK3 scheme. Fractional implicit Euler steps are once again used to integrate the Rayleigh-damping terms of the LHS sponge. This step is equivalent to the first one, except for the differences indicated in the compact description above.

 1. The right-hand sides are integrated over over $\Delta t / 2$ with an implicit Euler step, followed by the Poisson equation being solved and a correction step being performed. The Rayleigh-damping terms are doubled, since they were left out in the explicit Euler step. This step is equivalent to the second one, except for the differences indicated in the compact description above.

## MSGWaM

At the beginning of each time step, the saturation scheme is applied via integration of the respective term in the phase-space wave-action density equation with an explicit Euler step. i.e.

$$\mathcal{N}^{n + 1} = \mathcal{N}^n + \Delta t \mathcal{S}_s^n.$$

The ray equations are then integrated with the low-storage RK3 scheme, following

$$\begin{align*}
    q^{x, m + 1} & = \Delta t \frac{\partial \Omega^{m, n}}{\partial k} + \alpha_\mathrm{RK}^m q^{x, m}, & x^{m + 1} & = x^m + \beta_\mathrm{RK}^m q^{x, m + 1},\\
    q^{y, m + 1} & = \Delta t \frac{\partial \Omega^{m, n}}{\partial l}  + \alpha_\mathrm{RK}^m q^{y, m}, & y^{m + 1} & = y^m + \beta_\mathrm{RK}^m q^{y, m + 1},\\
    q^{z, m + 1} & = \Delta t \frac{\partial \Omega^{m, n}}{\partial m} + \alpha_\mathrm{RK}^m q^{z, m}, & z^{m + 1} & = z^m + \beta_\mathrm{RK}^m q^{z, m + 1},\\
    q^{k, m + 1} & = \Delta t \left(\frac{\partial \Omega^{m, n}}{\partial \widehat{x}} + G^{13} \frac{\partial \Omega^{m, n}}{\partial \widehat{z}}\right) + \alpha_\mathrm{RK}^m q^{k, m}, & k^{m + 1} & = k^m + \beta_\mathrm{RK}^m q^{k, m + 1},\\
    q^{l, m + 1} & = \Delta t \left(\frac{\partial \Omega^{m, n}}{\partial \widehat{y}} + G^{23} \frac{\partial \Omega^{m, n}}{\partial \widehat{z}}\right) + \alpha_\mathrm{RK}^m q^{l, m}, & l^{m + 1} & = l^m + \beta_\mathrm{RK}^m q^{l, m + 1},\\
    q^{m, m + 1} & = \frac{\Delta t}{J} \frac{\partial \Omega^{m, n}}{\partial \widehat{z}} + \alpha_\mathrm{RK}^m q^{m, m}, & m^{m + 1} & = m^m + \beta_\mathrm{RK}^m q^{m, m + 1},
\end{align*}$$

where $\Omega^{m, n} = \Omega \left(\boldsymbol{x}^m, \boldsymbol{k}^m, \widehat{\boldsymbol{u}}^n\right)$. At the end of every RK3 stage, the Rayleigh damping of the LHS sponge is applied to the phase-space wave-action density with a fractional implicit Euler step, i.e.

$$\mathcal{N}^{n + 1} \rightarrow \left(1 + 2 \alpha_\mathrm{R}^{n + 1} f_\mathrm{RK}^m \Delta t\right)^{- 1} \mathcal{N}^{n + 1}.$$

Note that the intrinsic frequency does not need to be updated since it is completely determined the other variables. Finally, the mean-flow tendencies for the current time step are calculated from the updated wave-property fields ([Jochum et al., 2025](https://doi.org/10.1175/JAS-D-24-0158.1)).

## Spatial discretization

PinCFlow is a finite-volume code that operates on a three-dimensional staggered C-grid ([Arakawa & Lamb, 1977](https://doi.org/10.1016/b978-0-12-460817-7.50009-4)), with scalar fields defined at the cell centers and velocity components defined at the respective interfaces ([Rieper et al., 2013](https://doi.org/10.1175/mwr-d-12-00026.1)). Advective-flux divergences are discretized with a monotone upwind scheme for conservation laws ([van Leer, 2003](https://doi.org/10.2514/6.2003-3559)) and a monotonized-centered variant limiter (e.g. [Kemm, 2010](https://doi.org/10.1002/fld.2357)). More specifically, the implementation is such that $P \widehat{\boldsymbol{u}}$ is the carrier flux (based on [Klein, 2009](https://doi.org/10.1007/s00162-009-0104-y); [Benacchio et al., 2014](https://doi.org/10.1175/mwr-d-13-00384.1); [Smolarkiewicz et al., 2014](https://doi.org/10.1016/j.jcp.2014.01.031) and [Benacchio & Klein, 2019](https://doi.org/10.1175/mwr-d-19-0073.1)). All other terms are discretized with centered differences (e.g. [Durran, 2010](https://doi.org/10.1007/978-1-4419-6412-0)). A Cartesian version of this is described in [Rieper et al. (2013)](https://doi.org/10.1175/mwr-d-12-00026.1) and [Schmid et al. (2021)](https://doi.org/10.1175/MWR-D-21-0126.1).

The spatial discretization of MSGWaM is based on the definition of ray volumes, six-dimensional cubes in phase space. The ray equations are integrated at the center points and surface midpoints, using trilinear interpolation to get mean-flow information at the locations of interest. To prevent uninhibited growing, ray volumes are split when their extent in any dimension of physical space exceeds the corresponding grid spacing. This is counteracted by merging in spectral space to keep the number of ray volumes below a user-defined threshold. In the computation of the mean-flow impact, the contribution of each ray volume is weighted by the fraction of the local grid cell covered by it. A more detailed description of the algorithm can be found in [Muraschko et al. (2014)](https://doi.org/10.1002/qj.2381), [Boeloeni et al. (2016)](https://doi.org/10.1175/JAS-D-16-0069.1), [Wilhelm et al. (2018)](https://doi.org/10.1175/JAS-D-17-0289.1), [Wei et al. (2019)](https://doi.org/10.1175/JAS-D-18-0337.1) and [Jochum et al. (2025)](https://doi.org/10.1175/JAS-D-24-0158.1).
