# Physics

## Coordinate system

PinCFlow uses a height-based terrain-following coordinate system ([Gal-Chen & Somerville, 1975](https://doi.org/10.1016/0021-9991(75)90037-6)) with vertical stretching. The transformation from the Cartesian system $\boldsymbol{x} = \left(x, y, z\right)^\mathrm{T}$ to the model's system $\widehat{\boldsymbol{x}} = \left(\widehat{x}, \widehat{y}, \widehat{z}\right)^\mathrm{T}$ is given by

$$\begin{align*}
    \widehat{z} & = \widehat{z} \left(\widetilde{z}\right), & \left(\widehat{x}, \widehat{y}, \widetilde{z}\right) & = \left(x, y, L_z \frac{z - h}{L_z - h}\right), &  \left(x, y, z\right) & = \left(\widehat{x}, \widehat{y}, \frac{L_z - h}{L_z}{\widetilde{z}} + h\right),
\end{align*}$$

where $L_z$ is the vertical extent of the model domain and $h$ is the surface topography. The contravariant basis vectors of the transformed system are

$$\begin{align*}
    \boldsymbol{\epsilon}^{\widehat{x}} & = \boldsymbol{e}^x, & \boldsymbol{\epsilon}^{\widehat{y}} & = \boldsymbol{e}^y, & \boldsymbol{\epsilon}^{\widehat{z}} = \frac{\partial h}{\partial \widehat{x}} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}} \boldsymbol{e}^x + \frac{\partial h}{\partial y} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}} \boldsymbol{e}^y + \frac{L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}}\boldsymbol{e}^z.
\end{align*}$$

The scalar products of these yield the metric tensor

$$G^{\mu \nu} = \begin{pmatrix}
    1 & 0 & \frac{\partial h}{\partial x} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}}\\
    0 & 1 & \frac{\partial h}{\partial y} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}}\\
    \frac{\partial h}{\partial x} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}} & \frac{\partial h}{\partial y} \frac{\widetilde{z} - L_z}{L_z - h} \frac{\partial \widehat{z}}{\partial \widetilde{z}} & \left\{\left(\frac{L_z}{L_z - h}\right)^2 + \left(\frac{\widetilde{z} - L_z}{L_z - h}\right)^2 \left[\left(\frac{\partial h}{\partial x}\right)^2 + \left(\frac{\partial h}{\partial y}\right)^2\right]\right\} \left(\frac{\partial \widehat{z}}{\partial \widetilde{z}}\right)^2
\end{pmatrix},$$

from which the Jacobian

$$J = \frac{L_z - h} {L_z} \frac{\partial \widetilde{z}}{\partial \widehat{z}}$$

can be determined. Finally, the transformation rule for the wind reads

$$\left(\widehat{u}, \widehat{v}, \widehat{w}\right) = \left(u, v, G^{13} u + G^{23} v + J^{- 1} w\right).$$

## Pseudo-incompressible equations

The atmosphere is decomposed into a hydrostatic background and deviations from it, specifically

$$\begin{align*}
    \pi \left(\boldsymbol{x}, t\right) & = \overline{\pi} \left(z\right) + \pi' \left(\boldsymbol{x}, t\right),\\
    \theta \left(\boldsymbol{x}, t\right) & = \overline{\theta} \left(z\right) + \theta' \left(\boldsymbol{x}, t\right),
\end{align*}$$

where $\pi$, $\theta$ and $t$ denote the Exner-pressure, potential temperature and time, respectively, and the background satisfies

$$\begin{align*}
    \overline{\pi} \left(0\right) = 1, && \frac{\mathrm{d} \overline{\pi}}{\mathrm{d} z} = - \frac{g}{c_p} \overline{\theta}^{- 1},
\end{align*}$$

with $g$ and $c_p$ being the gravitational acceleration and specific heat capacity at constant pressure, respectively ([Benacchio & Klein, 2019](https://doi.org/10.1175/mwr-d-19-0073.1)). Using this decomposition, the model equations in pseudo-incompressible mode are written as

$${\small\begin{align*}
    \frac{\partial \rho}{\partial t} + \frac{1}{J} \left(\frac{\partial J \rho u}{\partial \widehat{x}} + \frac{\partial J \rho v}{\partial \widehat{y}} + \frac{\partial J \rho \widehat{w}}{\partial \widehat{z}}\right) + \alpha_\mathrm{R} \left(\rho - \overline{\rho}\right) & = 0,\\
    \frac{\partial \rho'}{\partial t} + \frac{1}{J} \left(\frac{\partial J \rho' u}{\partial \widehat{x}} + \frac{\partial J \rho' v}{\partial \widehat{y}} + \frac{\partial J \rho' \widehat{w}}{\partial \widehat{z}}\right) + \alpha_\mathrm{R} \rho' & = \frac{N^2 \overline{\rho} w}{g},\\
    \frac{1}{J} \left(\frac{\partial J P u}{\partial \widehat{x}} + \frac{\partial J P v}{\partial \widehat{y}} + \frac{\partial J P \widehat{w}}{\partial \widehat{z}}\right) & = 0,\\
    \frac{\partial \rho u}{\partial t} + \mathcal{A}^{\rho u} - \mathcal{V}^{\rho u} - \mathcal{X}^{\rho u} - f \rho v + \alpha_\mathrm{R} \left(u - u_\mathrm{R}\right) & = - c_p P \mathcal{P}^{\rho u} - \beta_\mathrm{R} \rho u + F^{\rho u},\\
    \frac{\partial \rho v}{\partial t} + \mathcal{A}^{\rho v} - \mathcal{V}^{\rho v} - \mathcal{X}^{\rho v} + f \rho u + \alpha_\mathrm{R} \left(v - v_\mathrm{R}\right) & = - c_p P \mathcal{P}^{\rho v} - \beta_\mathrm{R} \rho v + F^{\rho v},\\
    \frac{\partial \rho \widehat{w}}{\partial t} + \mathcal{A}^{\rho \widehat{w}} - \mathcal{V}^{\rho \widehat{w}} - \mathcal{X}^{\rho \widehat{w}} - G^{13} f \rho v + G^{23} f \rho u + \alpha_\mathrm{R} \left(\widehat{w} - \widehat{w}_\mathrm{R}\right) & = - c_p P \mathcal{P}^{\rho \widehat{w}} - \frac{g \rho'}{J} - \beta_\mathrm{R} \rho \widehat{w} + F^{\rho \widehat{w}},
\end{align*}}$$

where $\rho \left(\boldsymbol{x}, t\right) = \overline{\rho} \left(z\right) + \rho' \left(\boldsymbol{x}, t\right)$ is the density, $P = \rho \theta = \overline{\rho} \overline{\theta}$ is the mass-weighted potential temperature, $N^2 = \left(g / \overline{\theta}\right) \left(\mathrm{d} \overline{\theta} / \mathrm{d} z\right)$ is the squared buoyancy frequency and $f = f_0$ is the Coriolis frequency. On the left-hand sides, $\alpha_\mathrm{R}$ and $\left(u_\mathrm{R}, v_\mathrm{R}, \widehat{w}_\mathrm{R}\right)^\mathrm{R}$ represent the Rayleigh-damping coefficient of a customizable sponge and the transformed wind that is to be obtained via the relaxation in it. In contrast, the Rayleigh-damping coefficient $\beta_\mathrm{R}$ on the right-hand side implements a preset sponge. The terms $\left(F^{\rho u}, F^{\rho v}, F^{\rho \widehat{w}}\right)^\mathrm{T}$ represent volume forces in the momentum equation, e.g. drag imposed by unresolved gravity waves. The advective momentum-flux divergences are given by

$$\begin{align*}
    \mathcal{A}^{\rho u} & = \frac{1}{J} \left(\frac{\partial J \rho u u}{\partial \widehat{x}} + \frac{\partial J \rho u v}{\partial \widehat{y}} + \frac{\partial J \rho u \widehat{w}}{\partial \widehat{z}}\right),\\
    \mathcal{A}^{\rho v} & = \frac{1}{J} \left(\frac{\partial J \rho v u}{\partial \widehat{x}} + \frac{\partial J \rho v v}{\partial \widehat{y}} + \frac{\partial J \rho v \widehat{w}}{\partial \widehat{z}}\right),\\
    \mathcal{A}^{\rho \widehat{w}} & = G^{13} \mathcal{A}^{\rho u} + G^{23} \mathcal{A}^{\rho v} + \frac{1}{J^2} \left(\frac{\partial J \rho w u}{\partial \widehat{x}} + \frac{\partial J \rho w v}{\partial \widehat{y}} + \frac{\partial J \rho w \widehat{w}}{\partial \widehat{z}}\right)
\end{align*}$$

and the viscous-flux divergences are

$$\begin{align*}
    \mathcal{V}^{\rho u} & = \frac{1}{J} \left(\frac{\partial J \eta \widehat{\Pi}^{11}}{\partial \widehat{x}} + \frac{\partial J \eta \widehat{\Pi}^{12}}{\partial \widehat{y}} + \frac{\partial J \eta \widehat{\Pi}^{13}}{\partial \widehat{z}}\right),\\
    \mathcal{V}^{\rho v} & = \frac{1}{J} \left(\frac{\partial J \eta \widehat{\Pi}^{21}}{\partial \widehat{x}} + \frac{\partial J \eta \widehat{\Pi}^{22}}{\partial \widehat{y}} + \frac{\partial J \eta \widehat{\Pi}^{23}}{\partial \widehat{z}}\right),\\
    \mathcal{V}^{\rho \widehat{w}} & = G^{13} \mathcal{V}^{\rho u} + G^{23} \mathcal{V}^{\rho v} + \frac{1}{J^2} \left[\frac{\partial J \eta \Pi^{31}}{\partial \widehat{x}} + \frac{\partial J \eta \Pi^{32}}{\partial \widehat{y}} + \frac{\partial}{\partial \widehat{z}} \left(J G^{13} \eta \Pi^{13} + J G^{23} \eta \Pi^{23} + \eta \Pi^{33}\right)\right].
\end{align*}$$

Therein, $\eta$ is the dynamic shear viscosity, the elements of the (symmetric) transformed stress tensor are given by

$$\begin{align*}
    \widehat{\Pi}^{11} & = \Pi^{11},\\
    \widehat{\Pi}^{12} & = \Pi^{12},\\
    \widehat{\Pi}^{22} & = \Pi^{22},\\
    \widehat{\Pi}^{13} & = G^{13} \Pi^{11} + G^{23} \Pi^{12} + \frac{1}{J} \Pi^{13},\\
    \widehat{\Pi}^{23} & = G^{13} \Pi^{12} + G^{23} \Pi^{22} + \frac{1}{J} \Pi^{23},\\
    \widehat{\Pi}^{33} & = \left(G^{13}\right)^2 \Pi^{11} + \left(G^{23}\right)^2 \Pi^{22} + \frac{1}{J^2} \Pi^{33} + 2 \left(G^{13} G^{23} \Pi^{12} + \frac{G^{13}}{J} \Pi^{13} + \frac{G^{23}}{J} \Pi^{23}\right)
\end{align*}$$

and those of its Cartesian counterpart are

$$\begin{align*}
    \Pi^{11} & = 2 \left(\frac{\partial u}{\partial \widehat{x}} + G^{13} \frac{\partial u}{\partial \widehat{z}}\right) - \frac{2}{3 J} \left(\frac{\partial J u}{\partial \widehat{x}} + \frac{\partial J v}{\partial \widehat{y}} + \frac{\partial J \widehat{w}}{\partial \widehat{z}}\right),\\
    \Pi^{12} & = \frac{\partial u}{\partial \widehat{y}} + G^{23} \frac{\partial u}{\partial \widehat{z}} + \frac{\partial v}{\partial \widehat{x}} + G^{13} \frac{\partial v}{\partial \widehat{z}},\\
    \Pi^{13} & = \frac{1}{J} \frac{\partial u}{\partial \widehat{z}} + \frac{\partial w}{{\partial \widehat{x}}} + G^{13} \frac{\partial w}{\partial \widehat{z}},\\
    \Pi^{22} & = 2 \left(\frac{\partial v}{{\partial \widehat{y}}} + G^{23} \frac{\partial v}{\partial \widehat{z}}\right) - \frac{2}{3 J} \left(\frac{\partial J u}{\partial \widehat{x}} + \frac{\partial J v}{\partial \widehat{y}} + \frac{\partial J \widehat{w}}{\partial \widehat{z}}\right),\\
    \Pi^{23} & = \frac{1}{J} \frac{\partial v}{\partial \widehat{z}} + \frac{\partial w}{\partial \widehat{y}} + G^{23} \frac{\partial w}{\partial \widehat{z}},\\
    \Pi^{33} & = \frac{2}{J} \frac{\partial w}{\partial \widehat{z}} - \frac{2}{3 J} \left(\frac{\partial J u}{\partial \widehat{x}} + \frac{\partial J v}{\partial \widehat{y}} + \frac{\partial J \widehat{w}}{\partial \widehat{z}}\right).
\end{align*}$$

In addition, artificial diffusion is represented by

$$\begin{align*}
    \mathcal{X}^{\rho u} & = \frac{1}{J} \left[\frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} u\right)}^{\widehat{x}}}{\partial \widehat{x}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} u\right)}^{\widehat{y}}}{\partial \widehat{y}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} u\right)}^{\widehat{z}}}{\partial \widehat{z}}\right],\\
    \mathcal{X}^{\rho v} & = \frac{1}{J} \left[\frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} v\right)}^{\widehat{x}}}{\partial \widehat{x}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} v\right)}^{\widehat{y}}}{\partial \widehat{y}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} v\right)}^{\widehat{z}}}{\partial \widehat{z}}\right),\\
    \mathcal{X}^{\rho \widehat{w}} & = G^{13} \mathcal{X}^{\rho u} + G^{23} \mathcal{X}^{\rho v} + \frac{1}{J^2} \left[\frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} w\right)}^{\widehat{x}}}{\partial \widehat{x}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} w\right)}^{\widehat{y}}}{\partial \widehat{y}} + \frac{\partial J \mu \widehat{\left(\boldsymbol{\nabla} w\right)}^{\widehat{z}}}{\partial \widehat{z}}\right],
\end{align*}$$

where $\mu$ is the momentum-diffusion coefficient and

$$\begin{align*}
    \widehat{\left(\boldsymbol{\nabla} \phi\right)}^{\widehat{x}} & = \frac{\partial \phi}{\partial \widehat{x}} + G^{13} \frac{\partial \phi}{\partial \widehat{z}},\\
    \widehat{\left(\boldsymbol{\nabla} \phi\right)}^{\widehat{y}} & = \frac{\partial \phi}{\partial \widehat{y}} + G^{23} \frac{\partial \phi}{\partial \widehat{z}},\\
    \widehat{\left(\boldsymbol{\nabla} \phi\right)}^{\widehat{z}} & = G^{13}\frac{\partial \phi}{\partial \widehat{x}} + G^{23} \frac{\partial \phi}{\partial \widehat{y}} + G^{33}\frac{\partial \phi}{\partial \widehat{z}}.
\end{align*}$$

Analogously, the components of the pressure gradient are given by

$$\begin{align*}
    \mathcal{P}^{\rho u} & = \frac{\partial \pi'}{\partial \widehat{x}} + G^{13} \frac{\partial \pi'}{\partial \widehat{z}},\\
    \mathcal{P}^{\rho v} & = \frac{\partial \pi'}{\partial \widehat{y}} + G^{23} \frac{\partial \pi'}{\partial \widehat{z}},\\
    \mathcal{P}^{\rho \widehat{w}} & = G^{13} \frac{\partial \pi'}{\partial \widehat{x}} + G^{23} \frac{\partial \pi'}{\partial \widehat{y}} + G^{33} \frac{\partial \pi'}{\partial \widehat{z}}
\end{align*}$$

 (see [Rieper et al., 2013](https://doi.org/10.1175/mwr-d-12-00026.1); [Schmid et al., 2021](https://doi.org/10.1175/MWR-D-21-0126.1)).

## Boussinesq equations

In Boussinesq mode, the continuity equation is removed and the density fluctuations are set to zero everywhere except in the auxiliary equation and the buoyancy term of the transformed-vertical-momentum equation. Moreover, $\overline{\rho}$, $\overline{\theta}$, $P$ and $N^2$ are replaced with the constant reference values $\rho_0$, $\theta_0$, $P_0$ and $N_0^2$.

## Compressible equations

In compressible mode, the identity $P = \overline{\rho} \overline{\theta}$ no longer holds, i.e. the mass-weighted potential temperature has a spatiotemporal dependence. The divergence constraint is thus replaced with the prognostic equation

$$\frac{\partial P}{\partial t} + \frac{1}{J} \left(\frac{\partial J P u}{\partial \widehat{x}} + \frac{\partial J P v}{\partial \widehat{y}} + \frac{\partial J P \widehat{w}}{\partial \widehat{z}}\right) - F^P + \alpha_\mathrm{R} P \left(1 - \frac{\overline{\rho}}{\rho}\right) = 0,$$

where the volume force $F^P$ represents a diabatic heating (e.g. due to heat conduction or unresolved gravity waves) that is not allowed in pseudo-incompressible mode. This term must also be represented in the auxiliary equation, which now reads

$$\frac{\partial \rho'}{\partial t} + \frac{1}{J} \left(\frac{\partial J \rho' u}{\partial \widehat{x}} + \frac{\partial J \rho' v}{\partial \widehat{y}} + \frac{\partial J \rho' \widehat{w}}{\partial \widehat{z}}\right) + \frac{F^P}{\overline{\theta}} + \alpha_\mathrm{R} \left[\rho' - \overline{\rho} \left(1 - \frac{P}{\rho \overline{\theta}}\right)\right] = \frac{N^2 P w}{g \overline{\theta}}.$$

Note that in addition to the new volume-force term, $\overline{\rho}$ has been replaced with $P / \overline{\theta}$, which is due to the density fluctuations being defined as $\rho' = \rho - P / \overline{\theta}$ ([Benacchio & Klein, 2019](https://doi.org/10.1175/mwr-d-19-0073.1)).

## Tracer transport

PinCFlow transports a passive tracer $\chi$ governed by

$$\frac{\partial \rho \chi}{\partial t} + \frac{1}{J} \left(\frac{\partial J \rho \chi u}{\partial \widehat{x}} + \frac{\partial J \rho \chi v}{\partial \widehat{y}} + \frac{\partial J \rho \chi \widehat{w}}{\partial \widehat{z}}\right) - F^{\rho \chi} + \alpha_\mathrm{R} \left[\rho \chi - \left(\rho \chi\right)^{\left(0\right)}\right] = 0,$$

where $F^{\rho \chi}$ represents another volume force, e.g. a tracer-flux divergence due to unresolved gravity waves, and $\left(\rho \chi\right)^{\left(0\right)}$ is the initial mass-weighted tracer.

## MSGWaM

### 3D transient theory

The gravity-wave dispersion relation is written as $\omega \left(\boldsymbol{x}, t\right) = \Omega \left[\boldsymbol{x}, t, \boldsymbol{k} \left(\boldsymbol{x}, t\right)\right]$, where in

$$\Omega \left(\boldsymbol{x}, t, \boldsymbol{k}\right) = \boldsymbol{k} \cdot \boldsymbol{u}_\mathrm{b} \pm \sqrt{\frac{N^2 \left(z\right) \left(k^2 + l^2\right) + f^2 m^2}{\left|\boldsymbol{k}\right|^2}},$$

the spatiotemporal dependence is split into an explicit part due to variations in the (resolved) background wind $\boldsymbol{u}_\mathrm{b} = \left(u_\mathrm{b}, v_\mathrm{b}, 0\right)^\mathrm{T}$ and the squared buoyancy frequency $N^2$, and an implicit part due to variations in the wavevector $\boldsymbol{k} = \left(k, l, m\right)^\mathrm{T}$. Using this decomposition, the eikonal equations are written as

$$\begin{align*}
    \dot{\omega} & = \left[\frac{\partial}{\partial t} + \boldsymbol{c}_\mathrm{g} \cdot \begin{pmatrix}
    \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
    \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
    J^{- 1} \partial_{\widehat{z}}
\end{pmatrix}\right] \omega = \frac{\partial \Omega}{\partial t},\\
\dot{\boldsymbol{k}} & = \left[\frac{\partial}{\partial t} + \boldsymbol{c}_\mathrm{g} \cdot \begin{pmatrix}
    \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
    \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
    J^{- 1} \partial_{\widehat{z}}
\end{pmatrix}\right] \boldsymbol{k} = \begin{pmatrix}
    \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
    \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
    J^{- 1} \partial_{\widehat{z}}
\end{pmatrix} \Omega.
\end{align*}$$

They are integrated along rays defined by the group velocity

$$\begin{align*}
    \dot{\boldsymbol{x}} = \boldsymbol{c}_\mathrm{g} = \boldsymbol{\nabla}_{\boldsymbol{k}} \Omega,
\end{align*}$$

where $\boldsymbol{\nabla}_{\boldsymbol{k}} = \left(\partial_k, \partial_l, \partial_m\right)^\mathrm{T}$. The gravity-wave energy is encoded in the phase-space wave-action density, which is governed by the equation

$$\dot{\mathcal{N}} = \left[\frac{\partial}{\partial t} + \boldsymbol{c}_\mathrm{g} \cdot \begin{pmatrix}
    \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
    \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
    J^{- 1} \partial_{\widehat{z}}
\end{pmatrix} + \dot{\boldsymbol{k}} \cdot \boldsymbol{\nabla}_{\boldsymbol{k}}\right] \mathcal{N} = \sum\limits_s \mathcal{S}_s,$$

where $\mathcal{S}_s$ are sinks and sources. This latter equation is integrated along rays defined by $\left(\dot{\boldsymbol{x}}, \dot{\boldsymbol{k}}\right)^\mathrm{T}$, so that $\mathcal{N}$ is conserved if the right-hand side is zero. The impact of the unresolved gravity waves on the resolved flow is given by

$$\begin{align*}
    \left(\frac{\partial \rho_\mathrm{b} u_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = - \frac{\rho_\mathrm{b}}{\overline{\rho}} \begin{pmatrix}
        \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
        \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
        J^{- 1} \partial_{\widehat{z}}
    \end{pmatrix} \cdot \left(\overline{\rho} \left\langle u' \boldsymbol{u}' \right\rangle\right) - \rho_\mathrm{b} \frac{f}{\overline{\theta}} \left\langle \theta' v' \right\rangle,\\
    \left(\frac{\partial \rho_\mathrm{b} v_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = - \frac{\rho_\mathrm{b}}{\overline{\rho}} \begin{pmatrix}
        \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
        \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
        J^{- 1} \partial_{\widehat{z}}
    \end{pmatrix} \cdot \left(\overline{\rho} \left\langle v' \boldsymbol{u}' \right\rangle\right) + \rho_\mathrm{b} \frac{f}{\overline{\theta}} \left\langle \theta' u' \right\rangle,\\
    \left(\frac{\partial \rho_\mathrm{b} \widehat{w}_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = G^{13} \left(\frac{\partial \rho_\mathrm{b} u_\mathrm{b}}{\partial t}\right)_\mathrm{w} + G^{23} \left(\frac{\partial \rho_\mathrm{b} v_\mathrm{b}}{\partial t}\right)_\mathrm{w},\\
    \left(\frac{\partial P_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = - \rho_\mathrm{b} \begin{pmatrix}
        \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
        \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
        0
    \end{pmatrix} \cdot \left\langle \theta' \boldsymbol{u}' \right\rangle,
\end{align*}$$

where $\rho_\mathrm{b}$, $\widehat{\boldsymbol{u}}_\mathrm{b} = \left(u_\mathrm{b}, v_\mathrm{b}, \widehat{w}_\mathrm{b}\right)^\mathrm{T}$ and $P_\mathrm{b}$ are the resolved density, transformed wind and mass-weighted potential temperature, respectively, and

$$\begin{align*}
    \overline{\rho} \left\langle u' u' \right\rangle & = \int \left[k \widehat{c}_{\mathrm{g} x} - \mathrm{sgn} \left(\left|f\right|\right) \frac{k \widehat{c}_{\mathrm{g} x} + l \widehat{c}_{\mathrm{g} y}}{1 - \left(\widehat{\omega} / f\right)^2}\right] \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \overline{\rho} \left\langle u' v' \right\rangle & = \int l \widehat{c}_{\mathrm{g} x} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \overline{\rho} \left\langle u' w' \right\rangle & = \int \frac{k \widehat{c}_{\mathrm{g} z}}{1 - \left(f / \widehat{\omega}\right)^2} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \overline{\rho} \left\langle v' v' \right\rangle & = \int \left[l \widehat{c}_{\mathrm{g} y} - \mathrm{sgn} \left(\left|f\right|\right) \frac{k \widehat{c}_{\mathrm{g} x} + l \widehat{c}_{\mathrm{g} y}}{1 - \left(\widehat{\omega} / f\right)^2}\right] \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \overline{\rho} \left\langle v' w' \right\rangle & = \int \frac{l \widehat{c}_{\mathrm{g} z}}{1 - \left(f / \widehat{\omega}\right)^2} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \left\langle \theta' u' \right\rangle & = \frac{f \overline{\theta}}{g \overline{\rho}} \int \frac{l m N^2}{\widehat{\omega} \left|\boldsymbol{k}\right|^2} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},\\
    \left\langle \theta' v' \right\rangle & = - \frac{f \overline{\theta}}{g \overline{\rho}} \int \frac{k m N^2}{\widehat{\omega} \left|\boldsymbol{k}\right|^2} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}},
\end{align*}$$

with $\widehat{\omega} = \omega - \boldsymbol{k} \cdot \boldsymbol{u}_\mathrm{b}$, $\widehat{\boldsymbol{c}}_\mathrm{g} = \left(\widehat{c}_{\mathrm{g} x}, \widehat{c}_{\mathrm{g} y}, \widehat{c}_{\mathrm{g} z}\right)^\mathrm{T} = \boldsymbol{\nabla}_{\boldsymbol{k}} \widehat{\omega}$ and $\mathrm{d} V_{\boldsymbol{k}} = \mathrm{d} k \mathrm{d} l \mathrm{d} m$ being the intrinsic frequency, intrinsic group velocity and spectral volume element, respectively (see [Achatz et al., 2017](https://doi.org/10.1002/qj.2926); [Achatz et al., 2023](https://doi.org/10.1063/5.0165180); [Jochum et al., 2025](https://doi.org/10.1175/JAS-D-24-0158.1)).

### Wave breaking

Wave breaking is captured with a saturation scheme. This scheme assumes that static instability leads to turbulent fluxes that may be described by a turbulent viscosity and diffusivity $K$ in a flux-gradient ansatz (see [Lindzen, 1981](https://doi.org/10.1029/JC086iC10p09707); [Becker, 2004](https://doi.org/10.1016/j.jastp.2004.01.019)). The divergence of these fluxes has a damping effect on the phase-space wave-action density, which is represented by the sink

$$\mathcal{S}_\mathrm{s} = - 2 K \left|\boldsymbol{k}\right|^2 \mathcal{N}.$$

This damping is assumed to be such that within one time step $\Delta t$, the instability criterion is no longer fulfilled, which implies

$$K = \frac{\overline{\rho}}{4 \Delta t} \left[\int N^4 \left(k^2 + l^2\right) m^2 \frac{\mathcal{N}}{\widehat{\omega}} \, \mathrm{d} V_{\boldsymbol{k}}\right]^{- 1} \max \left[0, \frac{2}{\overline{\rho}} \int \frac{N^4 \left(k^2 + l^2\right) m^2}{\widehat{\omega} \left|\boldsymbol{k}\right|^2} \mathcal{N} \mathrm{d} V_{\boldsymbol{k}} - \alpha_\mathrm{s}^2 N^4\right],$$

where $\alpha_\mathrm{s}$ is a saturation coefficient that accounts for uncertainties of the criterion ([Boeloeni et al., 2016](https://doi.org/10.1175/JAS-D-16-0069.1); [Boeloeni et al., 2021](https://doi.org/10.1175/JAS-D-20-0065.1)).

### Rayleigh damping

The Rayleigh damping in the LHS sponge introduced above is represented by the sink

$$\mathcal{S}_\mathrm{R} = - 2 \alpha_\mathrm{R} \mathcal{N}$$

in the phase-space-wave-action-density equation ([Jochum et al., 2025](https://doi.org/10.1175/JAS-D-24-0158.1)).

### 1D transient theory

The 3D transient theory can be reduced to a 1D transient one by removing all horizontal derivatives (in a Cartesian sense) and setting the horizontal components of the group velocity to zero. The eikonal equations then become

$$\begin{align*}
    \dot{\omega} & = \left(\frac{\partial}{\partial t} + \frac{c_{\mathrm{g} z}}{J} \frac{\partial}{\partial \widehat{z}}\right) \omega = \frac{\partial \Omega}{\partial t},\\
    \dot{m} & = \left(\frac{\partial}{\partial t} + \frac{c_{\mathrm{g} z}}{J} \frac{\partial}{\partial \widehat{z}}\right) m = J^{- 1} \frac{\partial \Omega}{\partial \widehat{z}}
\end{align*}$$

and are integrated along rays defined by the vertical group velocity

$$\dot{z} = c_{\mathrm{g} z} = \frac{\partial \Omega}{\partial m}.$$

The phase-space wave-action density equation is reduced to

$$\dot{\mathcal{N}} = \left(\frac{\partial}{\partial t} + \frac{c_{\mathrm{g} z}}{J} \frac{\partial}{\partial \widehat{z}} + \dot{m} \frac{\partial}{\partial m}\right) \mathcal{N} = \sum\limits_s \mathcal{S}_s$$

and integrated along rays defined by $\left(\dot{z}, \dot{m}\right)^\mathrm{T}$. Finally, the impact on the resolved flow becomes

$$\begin{align*}
    \left(\frac{\partial \rho_\mathrm{b} u_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = - \frac{\rho_\mathrm{b}}{J \overline{\rho}} \frac{\partial}{\partial \widehat{z}} \left(\overline{\rho} \left\langle u' w' \right\rangle\right) - \rho_\mathrm{b} \frac{f}{\overline{\theta}} \left\langle \theta' v' \right\rangle,\\
    \left(\frac{\partial \rho_\mathrm{b} v_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = - \frac{\rho_\mathrm{b}}{J \overline{\rho}} \frac{\partial}{\partial \widehat{z}} \left(\overline{\rho} \left\langle v' w' \right\rangle\right) + \rho_\mathrm{b} \frac{f}{\overline{\theta}} \left\langle \theta' u' \right\rangle,\\
    \left(\frac{\partial \rho_\mathrm{b} \widehat{w}_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = G^{13} \left(\frac{\partial \rho_\mathrm{b} u_\mathrm{b}}{\partial t}\right)_\mathrm{w} + G^{23} \left(\frac{\partial \rho_\mathrm{b} v_\mathrm{b}}{\partial t}\right)_\mathrm{w},\\
    \left(\frac{\partial P_\mathrm{b}}{\partial t}\right)_\mathrm{w} & = 0
\end{align*}$$

(e.g. [Boeloeni et al., 2016](https://doi.org/10.1175/JAS-D-16-0069.1); [Boeloeni et al., 2021](https://doi.org/10.1175/JAS-D-20-0065.1)).

### 1D steady-state theory

In the 1D steady-state theory, the phase-space-wave-action-density equation from the previous section is integrated over spectral space, the temporal derivative is removed and the quasilinear approximation is used to express the result in terms of the physical-space wave-action densities of individual spectral modes ($\mathcal{A}_\alpha$). This yields

$$J^{- 1} \frac{\partial}{\partial \widehat{z}} \left(c_{\mathrm{g} z, \alpha} \mathcal{A}_\alpha\right) = \sum\limits_s \mathcal{Q}_{s, \alpha},$$

where $\mathcal{N} = \sum_\alpha \mathcal{A}_\alpha \delta \left(\boldsymbol{k} - \boldsymbol{k}_\alpha\right)$ and $\mathcal{S}_s = \sum_\alpha \mathcal{Q}_{s, \alpha} \delta \left(\boldsymbol{k} - \boldsymbol{k}_\alpha\right)$. Similarly, the eikonal equations become

$$\begin{align*}
    \frac{c_{\mathrm{g} z, \alpha}}{J} \frac{\partial \omega_\alpha}{\partial \widehat{z}} & = 0,\\
    \frac{c_{\mathrm{g} z, \alpha}}{J} \frac{\partial m_\alpha}{\partial \widehat{z}} & = J^{- 1} \frac{\partial \Omega}{\partial \widehat{z}}.
\end{align*}$$

Thus, the intrinsic frequency and vertical wavenumber can be determined from

$$\begin{align*}
    \widehat{\omega}_\alpha & = \omega_\alpha - \boldsymbol{k}_\alpha \cdot \boldsymbol{u}_\mathrm{b},\\
    m_\alpha & = - \mathrm{sgn} \left(\widehat{\omega}_\alpha\right) \sqrt{\frac{\left(k_\alpha^2 + l_\alpha^2\right) \left(N^2 - \widehat{\omega}_\alpha^2\right)}{\widehat{\omega}_\alpha^2 - f^2}},
\end{align*}$$

where $\omega_\alpha$ is a constant. Since the wave-action-density equation is solved by vertical integration, the saturation scheme must be amended. Specifically, the saturation sink term $\mathcal{Q}_{\mathrm{s}, \alpha}$ is integrated over a pseudo-time step $J \Delta \widehat{z} / c_{\mathrm{g} z, \alpha}$. Therein, the turbulent viscosity and diffusivity is given by

$$K = \frac{\overline{\rho}}{4} \left[\sum\limits_\alpha \frac{J \Delta \widehat{z}}{c_{\mathrm{g} z, \alpha}} N^4 \left(k_\alpha^2 + l_\alpha^2\right) m_\alpha^2 \frac{\mathcal{A}_\alpha}{\widehat{\omega}_\alpha}\right]^{- 1} \max \left[0, \frac{2}{\overline{\rho}} \sum\limits_\alpha \frac{N^4 \left(k_\alpha^2 + l_\alpha^2\right) m_\alpha^2}{\widehat{\omega}_\alpha \left|\boldsymbol{k}_\alpha\right|^2} \mathcal{A}_\alpha - \alpha_\mathrm{s}^2 N^4\right]$$

(see [Boeloeni et al., 2021](https://doi.org/10.1175/JAS-D-20-0065.1); [Jochum et al., 2025](https://doi.org/10.1175/JAS-D-24-0158.1)).

### Orographic source

The orographic source is a lower-boundary condition for the wave-property fields. Its formulation uses the decomposition

$$h \left(x, y\right) = h_\mathrm{b} + \sum\limits_\alpha \Re \left\{h_{\mathrm{w}, \alpha} \exp \left[i \varphi_\alpha \left(x, y\right)\right]\right\}$$

of the surface topography $h$, where $h_\mathrm{b}$ and $h_\mathrm{w}$ vary slowly in $x$ and $y$, as opposed to the quickly varying orographic phase $\varphi_\alpha$. By inserting this and the WKB ansatz for the wind into the no-normal-flow condition

$$0 = \boldsymbol{u} \cdot \boldsymbol{n} = - u \frac{\partial h}{\partial x} - v \frac{\partial h}{\partial y} + w \qquad \mathrm{at} \quad z = h,$$

one obtains

$$\mathcal{N} = \frac{\overline{\rho}}{2} \sum\limits_\alpha \frac{\widehat{\omega}_\alpha \left|\boldsymbol{k}_\alpha\right|^2}{k_\alpha^2 + l_\alpha^2} \left|h_{\mathrm{w}, \alpha}\right|^2 \delta \left(\boldsymbol{k} - \boldsymbol{k}_\alpha\right) \qquad \mathrm{at} \quad z = h_\mathrm{b}.$$

The vertical wavenumber at the source is

$$m_\alpha = - \mathrm{sgn} \left(\widehat{\omega}_\alpha\right) \sqrt{\frac{\left(k_\alpha^2 + l_\alpha^2\right) \left(N^2 - \widehat{\omega}_\alpha^2\right)}{\widehat{\omega}_\alpha^2 - f^2}}$$

(see [Jochum et al., 2025](https://doi.org/10.1175/JAS-D-24-0158.1)).

### Tracer fluxes

The leading-order tracer fluxes due to unresolved gravity waves are given by

$$\begin{align*}
    \overline{\rho} \left\langle \chi' \boldsymbol{u}' \right\rangle & = f \int \frac{m}{\widehat{\omega} \left|\boldsymbol{k}\right|^2} \boldsymbol{k} \times \begin{pmatrix}
        \partial_{\widehat{x}} \chi_\mathrm{b} + G^{13} \partial_{\widehat{z}} \chi_\mathrm{b}\\
        \partial_{\widehat{y}} \chi_\mathrm{b} + G^{23} \partial_{\widehat{z}} \chi_\mathrm{b}\\
        J^{- 1} \partial_{\widehat{z}} \chi_\mathrm{b}
    \end{pmatrix} \mathcal{N} \, \mathrm{d} V_{\boldsymbol{k}}
\end{align*}$$

and the corresponding impact on the large-scale tracer is

$$\begin{align*}
    \left(\frac{\partial \rho_\mathrm{b} \chi_\mathrm{b}}{\partial t}\right)_\mathrm{w} = - \frac{\rho_\mathrm{b}}{\overline{\rho}} \begin{pmatrix}
        \partial_{\widehat{x}} + G^{13} \partial_{\widehat{z}}\\
        \partial_{\widehat{y}} + G^{23} \partial_{\widehat{z}}\\
        J^{- 1} \partial_{\widehat{z}}
    \end{pmatrix} \cdot \left(\overline{\rho} \left\langle \chi' \boldsymbol{u}' \right\rangle\right).
\end{align*}$$
