"""
```julia
compute_fluxes!(state::State, predictands::Predictands)
```

Compute fluxes by dispatching to specialized methods for each prognostic variable.

```julia
compute_fluxes!(state::State, predictands::Predictands, variable::Rho)
```

Compute the density fluxes in all three directions.

The fluxes are given by

```math
\\begin{align*}
    \\mathcal{F}^{\\rho, \\widehat{x}}_{i + 1 / 2} & = \\frac{\\tau_{\\widehat{x}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}^\\mathrm{R} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1}^\\mathrm{L}\\right\\},\\\\
    \\mathcal{F}^{\\rho, \\widehat{y}}_{j + 1 / 2} & = \\frac{\\tau_{\\widehat{y}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}^\\mathrm{F} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1}^\\mathrm{B}\\right\\},\\\\
    \\mathcal{F}^{\\rho, \\widehat{z}}_{k + 1 / 2} & = \\frac{\\tau_{\\widehat{z}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}^\\mathrm{U} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{k + 1}^\\mathrm{D}\\right\\},
\\end{align*}
```

where


```math
\\begin{align*}
    \\tau_{\\widehat{x}} & = \\left(J P_\\mathrm{old}\\right)_{i + 1 / 2} u_{\\mathrm{old}, i + 1 / 2},\\\\
    \\tau_{\\widehat{y}} & = \\left(J P_\\mathrm{old}\\right)_{j + 1 / 2} v_{\\mathrm{old}, j + 1 / 2},\\\\
    \\tau_{\\widehat{z}} & = \\left(J P_\\mathrm{old}\\right)_{k + 1 / 2} \\widehat{w}_{\\mathrm{old}, k + 1 / 2}
\\end{align*}
```

are the transporting velocities (weighted by the Jacobian) and ``\\widetilde{\\phi}`` is the reconstruction of ``\\rho / P_\\mathrm{old}``. More specifically, the superscripts ``\\mathrm{R}``, ``\\mathrm{L}``, ``\\mathrm{F}``, ``\\mathrm{B}``, ``\\mathrm{U}`` and ``\\mathrm{D}`` indicate reconstructions at the right, left, forward, backward, upward and downward cell interfaces of the respective grid points, respectively. Quantities with the subscript ``\\mathrm{old}`` are obtained from a previous state, which is partially passed to the method via `predictands`.


```julia
compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)
```

Compute the density-fluctuations fluxes in all three directions.

The computation is analogous to that of the density fluxes.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Union{Boussinesq, PseudoIncompressible},
    variable::P,
)
```

Return in non-compressible modes.

```julia
compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Compressible,
    variable::P,
)
```

Compute the mass-weighted potential-temperature fluxes in all three directions.

The fluxes are given by

```math
\\begin{align*}
    \\mathcal{F}^{P, \\widehat{x}}_{i + 1 / 2} & = \\left(J P_\\mathrm{old}\\right)_{i + 1 / 2} u_{\\mathrm{old}, i + 1 / 2},\\\\
    \\mathcal{F}^{P, \\widehat{y}}_{j + 1 / 2} & = \\left(J P_\\mathrm{old}\\right)_{j + 1 / 2} v_{\\mathrm{old}, j + 1 / 2},\\\\
    \\mathcal{F}^{P, \\widehat{z}}_{k + 1 / 2} & = \\left(J P_\\mathrm{old}\\right)_{k + 1 / 2} \\widehat{w}_{\\mathrm{old}, k + 1 / 2}.
\\end{align*}
```

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::U)
```

Compute the zonal-momentum fluxes in all three directions.

The fluxes are first set to the advective parts

```math
\\begin{align*}
    \\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} & = \\frac{\\tau_{\\widehat{x}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1 / 2}^\\mathrm{R} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{i + 3 / 2}^\\mathrm{L}\\right\\},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} & = \\frac{\\tau_{\\widehat{y}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1 / 2}^\\mathrm{F} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1 / 2, j + 1}^\\mathrm{B}\\right\\},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} & = \\frac{\\tau_{\\widehat{z}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1 / 2}^\\mathrm{U} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1 / 2, k + 1}^\\mathrm{D}\\right\\},
\\end{align*}
```

with

```math
\\begin{align*}
    \\tau_{\\widehat{x}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{i + 1 / 2} u_{\\mathrm{old}, i + 1 / 2}\\right]_{i + 1},\\\\
    \\tau_{\\widehat{y}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{j + 1 / 2} v_{\\mathrm{old}, j + 1 / 2}\\right]_{i + 1 / 2, j + 1 / 2},\\\\
    \\tau_{\\widehat{z}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{k + 1 / 2} \\widehat{w}_{\\mathrm{old}, k + 1 / 2}\\right]_{i + 1 / 2, k + 1 / 2}
\\end{align*}
```

and ``\\widetilde{\\phi}`` being the reconstruction of ``\\rho_{i + 1 / 2} u_{i + 1 / 2} / P_{\\mathrm{old}, i + 1 / 2}``. If the viscosity is nonzero, the viscous parts (weighted by the Jacobian) are then added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} - \\eta_{i + 1} \\left(J \\widehat{\\Pi}^{11}\\right)_{i + 1},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} - \\eta_{i + 1 / 2, j + 1 / 2} \\left(J \\widehat{\\Pi}^{12}\\right)_{i + 1 / 2, j + 1 / 2},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} - \\eta_{i + 1 / 2, k + 1 / 2} \\left(J \\widehat{\\Pi}^{13}\\right)_{i + 1 / 2, k + 1 / 2}.
\\end{align*}
```

Finally, if the diffusivity ``\\mu`` is nonzero, the diffusive parts (weighted by the Jacobian) are added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} - \\mu_{i + 1} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{x}}\\right]_{i + 1},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} - \\mu_{i + 1 / 2, j + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{y}}\\right]_{i + 1 / 2, j + 1 / 2},\\\\
    \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} - \\mu_{i + 1 / 2, k + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} u\\right)}^{\\widehat{z}}\\right]_{i + 1 / 2, k + 1 / 2}.
\\end{align*}
```

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::V)
```

Compute the meridional-momentum fluxes in all three directions.

The fluxes are first set to the advective parts

```math
\\begin{align*}
    \\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} & = \\frac{\\tau_{\\widehat{x}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1 / 2}^\\mathrm{R} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1, j + 1 / 2}^\\mathrm{L}\\right\\},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} & = \\frac{\\tau_{\\widehat{y}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1 / 2}^\\mathrm{F} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{j + 3 / 2}^\\mathrm{B}\\right\\},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} & = \\frac{\\tau_{\\widehat{z}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1 / 2}^\\mathrm{U} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1 / 2, k + 1}^\\mathrm{D}\\right\\},
\\end{align*}
```

with

```math
\\begin{align*}
    \\tau_{\\widehat{x}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{i + 1 / 2} u_{\\mathrm{old}, i + 1 / 2}\\right]_{i + 1 / 2, j + 1 / 2},\\\\
    \\tau_{\\widehat{y}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{j + 1 / 2} v_{\\mathrm{old}, j + 1 / 2}\\right]_{j + 1},\\\\
    \\tau_{\\widehat{z}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{k + 1 / 2} \\widehat{w}_{\\mathrm{old}, k + 1 / 2}\\right]_{j + 1 / 2, k + 1 / 2}
\\end{align*}
```

and ``\\widetilde{\\phi}`` being the reconstruction of ``\\rho_{j + 1 / 2} v_{j + 1 / 2} / P_{\\mathrm{old}, j + 1 / 2}``. If the viscosity is nonzero, the viscous parts (weighted by the Jacobian) are then added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} - \\eta_{i + 1 / 2, j + 1 / 2} \\left(J \\widehat{\\Pi}^{12}\\right)_{i + 1 / 2, j + 1 / 2},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} - \\eta_{j + 1} \\left(J \\widehat{\\Pi}^{22}\\right)_{j + 1},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} - \\eta_{j + 1 / 2, k + 1 / 2} \\left(J \\widehat{\\Pi}^{23}\\right)_{j + 1 / 2, k + 1 / 2}.
\\end{align*}
```

Finally, if the diffusivity ``\\mu`` is nonzero, the diffusive parts (weighted by the Jacobian) are added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} - \\mu_{i + 1 / 2, j + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{x}}\\right]_{i + 1 / 2, j + 1 / 2},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} - \\mu_{j + 1} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{y}}\\right]_{j + 1},\\\\
    \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} - \\mu_{j + 1 / 2, k + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} v\\right)}^{\\widehat{z}}\\right]_{j + 1 / 2, k + 1 / 2}.
\\end{align*}
```

```julia
compute_fluxes!(state::State, old_predictands::Predictands, variable::W)
```

Compute the vertical-momentum fluxes in all three directions.

The fluxes are first set to the advective parts

```math
\\begin{align*}
    \\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} & = \\frac{\\tau_{\\widehat{x}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{k + 1 / 2}^\\mathrm{R} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{x}}\\right)\\right] {\\widetilde{\\phi}}_{i + 1, k + 1 / 2}^\\mathrm{L}\\right\\},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} & = \\frac{\\tau_{\\widehat{y}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{k + 1 / 2}^\\mathrm{F} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{y}}\\right)\\right] {\\widetilde{\\phi}}_{j + 1, k + 1 / 2}^\\mathrm{B}\\right\\},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} & = \\frac{\\tau_{\\widehat{z}}}{2} \\left\\{\\left[1 + \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{k + 1 / 2}^\\mathrm{U} + \\left[1 - \\mathrm{sgn} \\left(\\tau_{\\widehat{z}}\\right)\\right] {\\widetilde{\\phi}}_{k + 3 / 2}^\\mathrm{D}\\right\\},
\\end{align*}
```

with

```math
\\begin{align*}
    \\tau_{\\widehat{x}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{i + 1 / 2} u_{\\mathrm{old}, i + 1 / 2}\\right]_{i + 1 / 2, k + 1 / 2},\\\\
    \\tau_{\\widehat{y}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{j + 1 / 2} v_{\\mathrm{old}, j + 1 / 2}\\right]_{j + 1 / 2, k + 1 / 2},\\\\
    \\tau_{\\widehat{z}} & = \\left[\\left(J P_\\mathrm{old}\\right)_{k + 1 / 2} \\widehat{w}_{\\mathrm{old}, k + 1 / 2}\\right]_{k + 1}
\\end{align*}
```

and ``\\widetilde{\\phi}`` being the reconstruction of ``\\rho_{k + 1 / 2} w_{k + 1 / 2} / P_{\\mathrm{old}, k + 1 / 2}``. If the viscosity is nonzero, the viscous parts (weighted by the Jacobian) are then added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} - \\eta_{i + 1 / 2, k + 1 / 2} \\left(J \\Pi^{13}\\right)_{i + 1 / 2, k + 1 / 2},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} - \\eta_{j + 1 / 2, k + 1 / 2} \\left(J \\Pi^{23}\\right)_{j + 1 / 2, k + 1 / 2},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} - \\eta_{k + 1} \\left[\\left(J G^{13} \\Pi^{13}\\right)_{k + 1} - \\left(J G^{23} \\Pi^{23}\\right)_{k + 1} - \\Pi^{33}_{k + 1}\\right].
\\end{align*}
```

Finally, if the diffusivity ``\\mu`` is nonzero, the diffusive parts (weighted by the Jacobian) are added, i.e.

```math
\\begin{align*}
    \\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} - \\mu_{i + 1 / 2, k + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{x}}\\right]_{i + 1 / 2, k + 1 / 2},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} - \\mu_{j + 1 / 2, k + 1 / 2} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{y}}\\right]_{j + 1 / 2, k + 1 / 2},\\\\
    \\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} & \\rightarrow \\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} - \\mu_{k + 1} \\left[J \\widehat{\\left(\\boldsymbol{\\nabla} w\\right)}^{\\widehat{z}}\\right]_{k + 1}.
\\end{align*}
```

```julia
compute_fluxes!(state::State, predictands::Predictands, tracer_setup::NoTracer)
```

Return for configurations without tracer transport.

```julia
compute_fluxes!(state::State, predictands::Predictands, tracer_setup::TracerOn)
```

Compute the tracer fluxes in all three directions.

The computation is analogous to that of the density fluxes.

```julia
compute_fluxes!(state::State, predictands::Predictands, variable::Theta)
```

Compute the potential temperature fluxes due to heat conduction (weighted by the Jacobian).

The fluxes are given by

```math
\\begin{align*}
    \\mathcal{F}^{\\theta, \\widehat{x}}_{i + 1 / 2} & = - \\lambda_{i + 1 / 2} \\left\\{\\frac{J_{i + 1 / 2}}{\\Delta \\widehat{x}} \\left[\\left(\\frac{P}{\\rho}\\right)_{i + 1} - \\frac{P}{\\rho}\\right]\\right.\\\\
    & \\qquad \\qquad \\qquad + \\left.\\frac{\\left(J G^{13}\\right)_{i + 1 / 2}}{2 \\Delta \\widehat{z}} \\left[\\left(\\frac{P}{\\rho}\\right)_{i + 1 / 2, k + 1} - \\left(\\frac{P}{\\rho}\\right)_{i + 1 / 2, k - 1}\\right]\\right\\},\\\\
    \\mathcal{F}^{\\theta, \\widehat{y}}_{j + 1 / 2} & = - \\lambda_{j + 1 / 2} \\left\\{\\frac{J_{j + 1 / 2}}{\\Delta \\widehat{y}} \\left[\\left(\\frac{P}{\\rho}\\right)_{j + 1} - \\frac{P}{\\rho}\\right]\\right.\\\\
    & \\qquad \\qquad \\qquad + \\left.\\frac{\\left(J G^{23}\\right)_{j + 1 / 2}}{2 \\Delta \\widehat{z}} \\left[\\left(\\frac{P}{\\rho}\\right)_{j + 1 / 2, k + 1} - \\left(\\frac{P}{\\rho}\\right)_{j + 1 / 2, k - 1}\\right]\\right\\},\\\\
    \\mathcal{F}^{\\theta, \\widehat{z}}_{k + 1 / 2} & = - \\lambda_{k + 1 / 2} \\left\\{\\frac{\\left(J G^{13}\\right)_{k + 1 / 2}}{2 \\Delta \\widehat{x}} \\left[\\left(\\frac{P}{\\rho}\\right)_{i + 1, k + 1 / 2} - \\left(\\frac{P}{\\rho}\\right)_{i - 1, k + 1 / 2}\\right]\\right.\\\\
    & \\qquad \\qquad \\qquad + \\frac{\\left(J G^{23}\\right)_{k + 1 / 2}}{2 \\Delta \\widehat{y}} \\left[\\left(\\frac{P}{\\rho}\\right)_{j + 1, k + 1 / 2} - \\left(\\frac{P}{\\rho}\\right)_{j - 1, k + 1 / 2}\\right]\\\\
    & \\qquad \\qquad \\qquad + \\left.\\frac{\\left(J G^{33}\\right)_{k + 1 / 2}}{\\Delta \\widehat{z}} \\left[\\left(\\frac{P}{\\rho}\\right)_{k + 1} - \\frac{P}{\\rho}\\right]\\right\\},
\\end{align*}
```

where ``\\lambda`` is the thermal conductivity (computed from `state.namelists.atmosphere.thermal_conductivity`).

# Arguments

  - `state`: Model state.

  - `predictands`/`old_predictands`: The predictands that are used to compute the transporting velocities.

  - `model`: Dynamic equations.

  - `variable`: Flux variable.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.FluxCalculator.compute_flux`](@ref)

  - [`PinCFlow.Update.compute_stress_tensor`](@ref)

  - [`PinCFlow.Update.conductive_heating`](@ref)
"""
function compute_fluxes! end

function compute_fluxes!(state::State, predictands::Predictands)
    (; model) = state.namelists.atmosphere

    compute_fluxes!(state, predictands, Rho())
    compute_fluxes!(state, predictands, RhoP())
    compute_fluxes!(state, predictands, U())
    compute_fluxes!(state, predictands, V())
    compute_fluxes!(state, predictands, W())

    compute_fluxes!(state, predictands, model, P())
    compute_fluxes!(state, predictands, state.namelists.tracer.tracer_setup)
    return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::Rho)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pbar, rhobar) = state.atmosphere
    (; rhotilde) = state.variables.reconstructions
    (; phirho) = state.variables.fluxes

    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhobaredger = 0.5 * (rhobar[i, j, k] + rhobar[i + 1, j, k])
        pedger = 0.5 * (pbar[i, j, k] + pbar[i + 1, j, k])
        rhor = rhotilde[i + 1, j, k, 1, 1] + rhobaredger / pedger
        rhol = rhotilde[i, j, k, 1, 2] + rhobaredger / pedger

        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        usurf = pedger * u0[i, j, k]

        frho = compute_flux(usurf, rhol, rhor)

        phirho[i, j, k, 1] = frho
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhobaredgef = 0.5 * (rhobar[i, j, k] + rhobar[i, j + 1, k])
        pedgef = 0.5 * (pbar[i, j, k] + pbar[i, j + 1, k])
        rhof = rhotilde[i, j + 1, k, 2, 1] + rhobaredgef / pedgef
        rhob = rhotilde[i, j, k, 2, 2] + rhobaredgef / pedgef

        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        vsurf = pedgef * v0[i, j, k]

        grho = compute_flux(vsurf, rhob, rhof)

        phirho[i, j, k, 2] = grho
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        rhobaredgeu =
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedgeu =
            (
                jac[i, j, k + 1] * pbar[i, j, k] +
                jac[i, j, k] * pbar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhou = rhotilde[i, j, k + 1, 3, 1] + rhobaredgeu / pedgeu
        rhod = rhotilde[i, j, k, 3, 2] + rhobaredgeu / pedgeu

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        wsurf = pedgeu * w0[i, j, k]

        hrho = compute_flux(wsurf, rhod, rhou)

        phirho[i, j, k, 3] = hrho
    end

    return
end

function compute_fluxes!(state::State, predictands::Predictands, variable::RhoP)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pbar) = state.atmosphere
    (; rhoptilde) = state.variables.reconstructions
    (; phirhop) = state.variables.fluxes

    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhor = rhoptilde[i + 1, j, k, 1, 1]
        rhol = rhoptilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        usurf = pedger * u0[i, j, k]

        frhop = compute_flux(usurf, rhol, rhor)

        phirhop[i, j, k, 1] = frhop
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhof = rhoptilde[i, j + 1, k, 2, 1]
        rhob = rhoptilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        vsurf = pedgef * v0[i, j, k]

        grhop = compute_flux(vsurf, rhob, rhof)

        phirhop[i, j, k, 2] = grhop
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        rhou = rhoptilde[i, j, k + 1, 3, 1]
        rhod = rhoptilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        wsurf = pedgeu * w0[i, j, k]

        hrhop = compute_flux(wsurf, rhod, rhou)

        phirhop[i, j, k, 3] = hrhop
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Union{Boussinesq, PseudoIncompressible},
    variable::P,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    model::Compressible,
    variable::P,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pbar) = state.atmosphere
    (; phip) = state.variables.fluxes

    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        phip[i, j, k, 1] =
            0.5 *
            (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            ) *
            u0[i, j, k]
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        phip[i, j, k, 2] =
            0.5 *
            (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            ) *
            v0[i, j, k]
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        phip[i, j, k, 3] =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1]) * w0[i, j, k]
    end

    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::U,
)
    (; grid) = state
    (; z_size) = state.namelists.domain
    (; re, uref, lref) = state.constants
    (; nz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pbar, rhobar) = state.atmosphere
    (; utilde) = state.variables.reconstructions
    (; phiu) = state.variables.fluxes
    (; kinematic_diffusivity) = state.namelists.atmosphere

    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    kmin = k0
    kmax = ko + nz == z_size ? k1 : k1 + 1

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in j0:j1, i in (i0 - 2):i1
        ur = utilde[i + 1, j, k, 1, 1]
        ul = utilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        predger =
            0.5 * (
                jac[i + 1, j, k] * pbar[i + 1, j, k] +
                jac[i + 2, j, k] * pbar[i + 2, j, k]
            )
        usurf = 0.5 * (pedger * u0[i, j, k] + predger * u0[i + 1, j, k])

        frhou = compute_flux(usurf, ul, ur)

        phiu[i, j, k, 1] = frhou
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        uf = utilde[i, j + 1, k, 2, 1]
        ub = utilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        predgef =
            0.5 * (
                jac[i + 1, j, k] * pbar[i + 1, j, k] +
                jac[i + 1, j + 1, k] * pbar[i + 1, j + 1, k]
            )
        vsurf = 0.5 * (pedgef * v0[i, j, k] + predgef * v0[i + 1, j, k])

        grhou = compute_flux(vsurf, ub, uf)

        phiu[i, j, k, 2] = grhou
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in j0:j1, i in (i0 - 1):i1
        uu = utilde[i, j, k + 1, 3, 1]
        ud = utilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        predgeu =
            jac[i + 1, j, k] *
            jac[i + 1, j, k + 1] *
            (pbar[i + 1, j, k] + pbar[i + 1, j, k + 1]) /
            (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + predgeu * w0[i + 1, j, k])

        hrhou = compute_flux(wsurf, ud, uu)

        phiu[i, j, k, 3] = hrhou
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps() && kinematic_diffusivity == 0.0
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in j0:j1, i in (i0 - 2):i1
        coef_v = 1 / re * rhobar[i + 1, j, k0]

        frhou_visc =
            coef_v *
            jac[i + 1, j, k] *
            compute_stress_tensor(i + 1, j, k, 1, 1, state)

        phiu[i, j, k, 1] -= frhou_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_v =
            1 / re *
            0.25 *
            (
                rhobar[i, j, k0] +
                rhobar[i + 1, j, k0] +
                rhobar[i, j + 1, k0] +
                rhobar[i + 1, j + 1, k0]
            )

        grhou_visc =
            coef_v *
            0.25 *
            (
                jac[i, j, k] * compute_stress_tensor(i, j, k, 1, 2, state) +
                jac[i + 1, j, k] *
                compute_stress_tensor(i + 1, j, k, 1, 2, state) +
                jac[i, j + 1, k] *
                compute_stress_tensor(i, j + 1, k, 1, 2, state) +
                jac[i + 1, j + 1, k] *
                compute_stress_tensor(i + 1, j + 1, k, 1, 2, state)
            )

        phiu[i, j, k, 2] -= grhou_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in j0:j1, i in (i0 - 1):i1
        coef_v = 1 / re * 0.5 * (rhobar[i, j, k0] + rhobar[i + 1, j, k0])

        stresstens13 =
            met[i, j, k, 1, 3] * compute_stress_tensor(i, j, k, 1, 1, state) +
            met[i, j, k, 2, 3] * compute_stress_tensor(i, j, k, 1, 2, state) +
            compute_stress_tensor(i, j, k, 1, 3, state) / jac[i, j, k]
        stresstens13r =
            met[i + 1, j, k, 1, 3] *
            compute_stress_tensor(i + 1, j, k, 1, 1, state) +
            met[i + 1, j, k, 2, 3] *
            compute_stress_tensor(i + 1, j, k, 1, 2, state) +
            compute_stress_tensor(i + 1, j, k, 1, 3, state) / jac[i + 1, j, k]
        stresstens13u =
            met[i, j, k + 1, 1, 3] *
            compute_stress_tensor(i, j, k + 1, 1, 1, state) +
            met[i, j, k + 1, 2, 3] *
            compute_stress_tensor(i, j, k + 1, 1, 2, state) +
            compute_stress_tensor(i, j, k + 1, 1, 3, state) / jac[i, j, k + 1]
        stresstens13ru =
            met[i + 1, j, k + 1, 1, 3] *
            compute_stress_tensor(i + 1, j, k + 1, 1, 1, state) +
            met[i + 1, j, k + 1, 2, 3] *
            compute_stress_tensor(i + 1, j, k + 1, 1, 2, state) +
            compute_stress_tensor(i + 1, j, k + 1, 1, 3, state) /
            jac[i + 1, j, k + 1]
        hrhou_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (stresstens13 + stresstens13u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (stresstens13r + stresstens13ru) /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiu[i, j, k, 3] -= hrhou_visc
    end

    #-------------------------------------------------------------------
    #             Diffusion fluxes
    #-------------------------------------------------------------------

    if kinematic_diffusivity == 0.0
        return
    end

    mu_mom_diff = kinematic_diffusivity / uref / lref

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in j0:j1, i in (i0 - 2):i1
        coef_d = mu_mom_diff * rhobar[i + 1, j, k0]

        frhou_diff =
            coef_d *
            jac[i + 1, j, k] *
            compute_momentum_diffusion_terms(state, i + 1, j, k, U(), X())

        phiu[i, j, k, 1] -= frhou_diff
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_d =
            mu_mom_diff *
            0.25 *
            (
                rhobar[i, j, k0] +
                rhobar[i + 1, j, k0] +
                rhobar[i, j + 1, k0] +
                rhobar[i + 1, j + 1, k0]
            )

        grhou_diff =
            coef_d *
            0.25 *
            (
                jac[i, j, k] *
                compute_momentum_diffusion_terms(state, i, j, k, U(), Y()) +
                jac[i + 1, j, k] *
                compute_momentum_diffusion_terms(state, i + 1, j, k, U(), Y()) +
                jac[i, j + 1, k] *
                compute_momentum_diffusion_terms(state, i, j + 1, k, U(), Y()) +
                jac[i + 1, j + 1, k] * compute_momentum_diffusion_terms(
                    state,
                    i + 1,
                    j + 1,
                    k,
                    U(),
                    Y(),
                )
            )

        phiu[i, j, k, 2] -= grhou_diff
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in j0:j1, i in (i0 - 1):i1
        coef_dr = mu_mom_diff * rhobar[i + 1, j, k0]

        coef_dl = mu_mom_diff * rhobar[i, j, k0]

        coef_d = 0.5 * (coef_dr + coef_dl)

        mom_diff =
            jac[i, j, k] *
            compute_momentum_diffusion_terms(state, i, j, k, U(), Z())

        mom_diff_r =
            jac[i + 1, j, k] *
            compute_momentum_diffusion_terms(state, i + 1, j, k, U(), Z())

        mom_diff_u =
            jac[i, j, k + 1] *
            compute_momentum_diffusion_terms(state, i, j, k + 1, U(), Z())

        mom_diff_ru =
            jac[i + 1, j, k + 1] *
            compute_momentum_diffusion_terms(state, i + 1, j, k + 1, U(), Z())

        hrhou_diff =
            coef_d *
            0.5 *
            (
                jac[i, j, k] * jac[i, j, k + 1] * (mom_diff + mom_diff_u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (mom_diff_r + mom_diff_ru) /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiu[i, j, k, 3] -= hrhou_diff
    end

    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::V,
)
    (; grid) = state
    (; z_size) = state.namelists.domain
    (; re, uref, lref) = state.constants
    (; nz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pbar, rhobar) = state.atmosphere
    (; vtilde) = state.variables.reconstructions
    (; phiv) = state.variables.fluxes
    (; kinematic_diffusivity) = state.namelists.atmosphere

    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    kmin = k0
    kmax = ko + nz == z_size ? k1 : k1 + 1

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        vr = vtilde[i + 1, j, k, 1, 1]
        vl = vtilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        pfedger =
            0.5 * (
                jac[i, j + 1, k] * pbar[i, j + 1, k] +
                jac[i + 1, j + 1, k] * pbar[i + 1, j + 1, k]
            )
        usurf = 0.5 * (pedger * u0[i, j, k] + pfedger * u0[i, j + 1, k])

        frhov = compute_flux(usurf, vl, vr)

        phiv[i, j, k, 1] = frhov
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 2):j1, i in i0:i1
        vf = vtilde[i, j + 1, k, 2, 1]
        vb = vtilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        pfedgef =
            0.5 * (
                jac[i, j + 1, k] * pbar[i, j + 1, k] +
                jac[i, j + 2, k] * pbar[i, j + 2, k]
            )
        vsurf = 0.5 * (pedgef * v0[i, j, k] + pfedgef * v0[i, j + 1, k])

        grhov = compute_flux(vsurf, vb, vf)

        phiv[i, j, k, 2] = grhov
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in (j0 - 1):j1, i in i0:i1
        vu = vtilde[i, j, k + 1, 3, 1]
        vd = vtilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        pfedgeu =
            jac[i, j + 1, k] *
            jac[i, j + 1, k + 1] *
            (pbar[i, j + 1, k] + pbar[i, j + 1, k + 1]) /
            (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + pfedgeu * w0[i, j + 1, k])

        hrhov = compute_flux(wsurf, vd, vu)

        phiv[i, j, k, 3] = hrhov
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps() && kinematic_diffusivity == 0.0
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_v =
            1 / re *
            0.25 *
            (
                rhobar[i, j, k0] +
                rhobar[i + 1, j, k0] +
                rhobar[i, j + 1, k0] +
                rhobar[i + 1, j + 1, k0]
            )

        frhov_visc =
            coef_v *
            0.25 *
            (
                jac[i, j, k] * compute_stress_tensor(i, j, k, 2, 1, state) +
                jac[i + 1, j, k] *
                compute_stress_tensor(i + 1, j, k, 2, 1, state) +
                jac[i, j + 1, k] *
                compute_stress_tensor(i, j + 1, k, 2, 1, state) +
                jac[i + 1, j + 1, k] *
                compute_stress_tensor(i + 1, j + 1, k, 2, 1, state)
            )

        phiv[i, j, k, 1] -= frhov_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 2):j1, i in i0:i1
        coef_v = 1 / re * rhobar[i, j + 1, k0]

        grhov_visc =
            coef_v *
            jac[i, j + 1, k] *
            compute_stress_tensor(i, j + 1, k, 2, 2, state)

        phiv[i, j, k, 2] -= grhov_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in (j0 - 1):j1, i in i0:i1
        coef_v = 1 / re * 0.5 * (rhobar[i, j, k0] + rhobar[i, j + 1, k0])

        stresstens23 =
            met[i, j, k, 1, 3] * compute_stress_tensor(i, j, k, 2, 1, state) +
            met[i, j, k, 2, 3] * compute_stress_tensor(i, j, k, 2, 2, state) +
            compute_stress_tensor(i, j, k, 2, 3, state) / jac[i, j, k]
        stresstens23f =
            met[i, j + 1, k, 1, 3] *
            compute_stress_tensor(i, j + 1, k, 2, 1, state) +
            met[i, j + 1, k, 2, 3] *
            compute_stress_tensor(i, j + 1, k, 2, 2, state) +
            compute_stress_tensor(i, j + 1, k, 2, 3, state) / jac[i, j + 1, k]
        stresstens23u =
            met[i, j, k + 1, 1, 3] *
            compute_stress_tensor(i, j, k + 1, 2, 1, state) +
            met[i, j, k + 1, 2, 3] *
            compute_stress_tensor(i, j, k + 1, 2, 2, state) +
            compute_stress_tensor(i, j, k + 1, 2, 3, state) / jac[i, j, k + 1]
        stresstens23fu =
            met[i, j + 1, k + 1, 1, 3] *
            compute_stress_tensor(i, j + 1, k + 1, 2, 1, state) +
            met[i, j + 1, k + 1, 2, 3] *
            compute_stress_tensor(i, j + 1, k + 1, 2, 2, state) +
            compute_stress_tensor(i, j + 1, k + 1, 2, 3, state) /
            jac[i, j + 1, k + 1]
        hrhov_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (stresstens23 + stresstens23u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (stresstens23f + stresstens23fu) /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiv[i, j, k, 3] -= hrhov_visc
    end

    #-------------------------------------------------------------------
    #                          Diffusion fluxes
    #-------------------------------------------------------------------

    if kinematic_diffusivity == 0.0
        return
    end

    mu_mom_diff = kinematic_diffusivity / uref / lref

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in (i0 - 1):i1
        coef_d =
            mu_mom_diff *
            0.25 *
            (
                rhobar[i, j, k0] +
                rhobar[i + 1, j, k0] +
                rhobar[i, j + 1, k0] +
                rhobar[i + 1, j + 1, k0]
            )

        frhov_diff =
            coef_d *
            0.25 *
            (
                jac[i, j, k] *
                compute_momentum_diffusion_terms(state, i, j, k, V(), X()) +
                jac[i + 1, j, k] *
                compute_momentum_diffusion_terms(state, i + 1, j, k, V(), X()) +
                jac[i, j + 1, k] *
                compute_momentum_diffusion_terms(state, i, j + 1, k, V(), X()) +
                jac[i + 1, j + 1, k] * compute_momentum_diffusion_terms(
                    state,
                    i + 1,
                    j + 1,
                    k,
                    V(),
                    X(),
                )
            )

        phiv[i, j, k, 1] -= frhov_diff
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in kmin:kmax, j in (j0 - 2):j1, i in i0:i1
        coef_d = mu_mom_diff * rhobar[i, j + 1, k0]

        grhov_diff =
            coef_d *
            jac[i, j + 1, k] *
            compute_momentum_diffusion_terms(state, i, j + 1, k, V(), Y())

        phiv[i, j, k, 2] -= grhov_diff
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (kmin - 1):kmax, j in (j0 - 1):j1, i in i0:i1
        coef_dr = mu_mom_diff * rhobar[i, j + 1, k0]

        coef_dl = mu_mom_diff * rhobar[i, j, k0]

        coef_d = 0.5 * (coef_dr + coef_dl)

        u_diff =
            jac[i, j, k] *
            compute_momentum_diffusion_terms(state, i, j, k, V(), Z())

        u_diff_f =
            jac[i, j + 1, k] *
            compute_momentum_diffusion_terms(state, i, j + 1, k, V(), Z())

        u_diff_u =
            jac[i, j, k + 1] *
            compute_momentum_diffusion_terms(state, i, j, k + 1, V(), Z())

        u_diff_fu =
            jac[i, j + 1, k + 1] *
            compute_momentum_diffusion_terms(state, i, j + 1, k + 1, V(), Z())

        hrhov_diff =
            coef_d *
            0.5 *
            (
                jac[i, j, k] * jac[i, j, k + 1] * (u_diff + u_diff_u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (u_diff_f + u_diff_fu) /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiv[i, j, k, 3] -= hrhov_diff
    end

    return
end

function compute_fluxes!(
    state::State,
    old_predictands::Predictands,
    variable::W,
)
    (; grid) = state
    (; re, uref, lref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = grid
    (; pbar, rhobar) = state.atmosphere
    (; wtilde) = state.variables.reconstructions
    (; phiw) = state.variables.fluxes
    (; kinematic_diffusivity) = state.namelists.atmosphere

    (u0, v0, w0) = (old_predictands.u, old_predictands.v, old_predictands.w)

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in (i0 - 1):i1
        wr = wtilde[i + 1, j, k, 1, 1]
        wl = wtilde[i, j, k, 1, 2]

        pedger =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        puedger =
            0.5 * (
                jac[i, j, k + 1] * pbar[i, j, k + 1] +
                jac[i + 1, j, k + 1] * pbar[i + 1, j, k + 1]
            )
        usurf =
            (
                (jac[i, j, k + 1] + jac[i + 1, j, k + 1]) *
                pedger *
                u0[i, j, k] +
                (jac[i, j, k] + jac[i + 1, j, k]) * puedger * u0[i, j, k + 1]
            ) / (
                jac[i, j, k] +
                jac[i + 1, j, k] +
                jac[i, j, k + 1] +
                jac[i + 1, j, k + 1]
            )

        frhow = compute_flux(usurf, wl, wr)

        phiw[i, j, k, 1] = frhow
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in (j0 - 1):j1, i in i0:i1
        wf = wtilde[i, j + 1, k, 2, 1]
        wb = wtilde[i, j, k, 2, 2]

        pedgef =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        puedgef =
            0.5 * (
                jac[i, j, k + 1] * pbar[i, j, k + 1] +
                jac[i, j + 1, k + 1] * pbar[i, j + 1, k + 1]
            )
        vsurf =
            (
                (jac[i, j, k + 1] + jac[i, j + 1, k + 1]) *
                pedgef *
                v0[i, j, k] +
                (jac[i, j, k] + jac[i, j + 1, k]) * puedgef * v0[i, j, k + 1]
            ) / (
                jac[i, j, k] +
                jac[i, j + 1, k] +
                jac[i, j, k + 1] +
                jac[i, j + 1, k + 1]
            )

        grhow = compute_flux(vsurf, wb, wf)

        phiw[i, j, k, 2] = grhow
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 2):k1, j in j0:j1, i in i0:i1
        wu = wtilde[i, j, k + 1, 3, 1]
        wd = wtilde[i, j, k, 3, 2]

        pedgeu =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        puedgeu =
            jac[i, j, k + 1] *
            jac[i, j, k + 2] *
            (pbar[i, j, k + 1] + pbar[i, j, k + 2]) /
            (jac[i, j, k + 1] + jac[i, j, k + 2])
        wsurf = 0.5 * (pedgeu * w0[i, j, k] + puedgeu * w0[i, j, k + 1])

        hrhow = compute_flux(wsurf, wd, wu)

        phiw[i, j, k, 3] = hrhow
    end

    #-------------------------------------------------------------------
    #                          Viscous fluxes
    #-------------------------------------------------------------------

    if 1 / re <= eps() && kinematic_diffusivity == 0.0
        return
    end

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in (i0 - 1):i1
        coef_v = 1 / re * 0.5 * (rhobar[i, j, k0] + rhobar[i + 1, j, k0])

        frhow_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (
                    compute_stress_tensor(i, j, k, 3, 1, state) +
                    compute_stress_tensor(i, j, k + 1, 3, 1, state)
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (
                    compute_stress_tensor(i + 1, j, k, 3, 1, state) +
                    compute_stress_tensor(i + 1, j, k + 1, 3, 1, state)
                ) / (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiw[i, j, k, 1] -= frhow_visc
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in (j0 - 1):j1, i in i0:i1
        coef_v = 1 / re * 0.5 * (rhobar[i, j, k0] + rhobar[i, j + 1, k0])

        grhow_visc =
            coef_v *
            0.5 *
            (
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (
                    compute_stress_tensor(i, j, k, 3, 1, state) +
                    compute_stress_tensor(i, j, k + 1, 3, 1, state)
                ) / (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (
                    compute_stress_tensor(i, j + 1, k, 3, 1, state) +
                    compute_stress_tensor(i, j + 1, k + 1, 3, 1, state)
                ) / (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiw[i, j, k, 2] -= grhow_visc
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 2):k1, j in j0:j1, i in i0:i1
        coef_v = 1 / re * rhobar[i, j, k0]

        hrhow_visc =
            coef_v * (
                jac[i, j, k + 1] *
                met[i, j, k + 1, 1, 3] *
                compute_stress_tensor(i, j, k + 1, 3, 1, state) +
                jac[i, j, k + 1] *
                met[i, j, k + 1, 2, 3] *
                compute_stress_tensor(i, j, k + 1, 3, 2, state) +
                compute_stress_tensor(i, j, k + 1, 3, 3, state)
            )

        phiw[i, j, k, 3] -= hrhow_visc
    end

    #-------------------------------------------------------------------
    #                          Diffusion fluxes
    #-------------------------------------------------------------------

    if kinematic_diffusivity == 0.0
        return
    end

    mu_mom_diff = kinematic_diffusivity / uref / lref

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in (i0 - 1):i1
        coef_dr = mu_mom_diff * rhobar[i + 1, j, k0]

        coef_dl = mu_mom_diff * rhobar[i, j, k0]

        coef_d = 0.5 * (coef_dr + coef_dl)

        w_diff =
            jac[i, j, k] *
            compute_momentum_diffusion_terms(state, i, j, k, W(), X())

        w_diff_r =
            jac[i + 1, j, k] *
            compute_momentum_diffusion_terms(state, i + 1, j, k, W(), X())

        w_diff_u =
            jac[i, j, k + 1] *
            compute_momentum_diffusion_terms(state, i, j, k + 1, W(), X())

        w_diff_ru =
            jac[i + 1, j, k + 1] *
            compute_momentum_diffusion_terms(state, i + 1, j, k + 1, W(), X())

        frhow_diff =
            coef_d *
            0.5 *
            (
                jac[i, j, k] * jac[i, j, k + 1] * (w_diff + w_diff_u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i + 1, j, k] *
                jac[i + 1, j, k + 1] *
                (w_diff_r + w_diff_ru) /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1])
            )

        phiw[i, j, k, 1] -= frhow_diff
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in (j0 - 1):j1, i in i0:i1
        coef_dr = mu_mom_diff * rhobar[i, j + 1, k0]

        coef_dl = mu_mom_diff * rhobar[i, j, k0]

        coef_d = 0.5 * (coef_dr + coef_dl)

        w_diff =
            jac[i, j, k] *
            compute_momentum_diffusion_terms(state, i, j, k, W(), Y())

        w_diff_f =
            jac[i, j + 1, k] *
            compute_momentum_diffusion_terms(state, i, j + 1, k, W(), Y())

        w_diff_u =
            jac[i, j, k + 1] *
            compute_momentum_diffusion_terms(state, i, j, k + 1, W(), Y())

        w_diff_fu =
            jac[i, j + 1, k + 1] *
            compute_momentum_diffusion_terms(state, i, j + 1, k + 1, W(), Y())

        grhow_diff =
            coef_d *
            0.5 *
            (
                jac[i, j, k] * jac[i, j, k + 1] * (w_diff + w_diff_u) /
                (jac[i, j, k] + jac[i, j, k + 1]) +
                jac[i, j + 1, k] *
                jac[i, j + 1, k + 1] *
                (w_diff_f + w_diff_fu) /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1])
            )

        phiw[i, j, k, 2] -= grhow_diff
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 2):k1, j in j0:j1, i in i0:i1
        coef_d = mu_mom_diff * rhobar[i, j, k0]

        hrhow_visc =
            coef_d *
            jac[i, j, k + 1] *
            compute_momentum_diffusion_terms(state, i, j, k + 1, W(), Z())

        phiw[i, j, k, 3] -= hrhow_visc
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    tracer_setup::NoTracer,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    tracer_setup::TracerOn,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pbar) = state.atmosphere
    (; tracerreconstructions, tracerfluxes) = state.tracer

    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    @ivy for field in 1:fieldcount(TracerPredictands)
        chir = getfield(tracerreconstructions, field)[2:end, :, :, 1, 1]
        chil = getfield(tracerreconstructions, field)[:, :, :, 1, 2]
        fchi = getfield(tracerfluxes, field)[:, :, :, 1]
        for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
            pedger =
                0.5 * (
                    jac[i, j, k] * pbar[i, j, k] +
                    jac[i + 1, j, k] * pbar[i + 1, j, k]
                )
            usurf = pedger * u0[i, j, k]

            fchi[i, j, k] = compute_flux(usurf, chil[i, j, k], chir[i, j, k])
        end

        chif = getfield(tracerreconstructions, field)[:, 2:end, :, 2, 1]
        chib = getfield(tracerreconstructions, field)[:, :, :, 2, 2]
        gchi = getfield(tracerfluxes, field)[:, :, :, 2]
        for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
            pedgef =
                0.5 * (
                    jac[i, j, k] * pbar[i, j, k] +
                    jac[i, j + 1, k] * pbar[i, j + 1, k]
                )
            vsurf = pedgef * v0[i, j, k]

            gchi[i, j, k] = compute_flux(vsurf, chib[i, j, k], chif[i, j, k])
        end

        chiu = getfield(tracerreconstructions, field)[:, :, 2:end, 3, 1]
        chid = getfield(tracerreconstructions, field)[:, :, :, 3, 2]
        hchi = getfield(tracerfluxes, field)[:, :, :, 3]
        for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
            pedgeu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (pbar[i, j, k] + pbar[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
            wsurf = pedgeu * w0[i, j, k]

            hchi[i, j, k] = compute_flux(wsurf, chid[i, j, k], chiu[i, j, k])
        end
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    variable::Theta,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, dx, dy, dz, met) = state.grid
    (; pbar, rhobar) = state.atmosphere
    (; phitheta) = state.variables.fluxes
    (; thermal_conductivity) = state.namelists.atmosphere
    (; uref, lref) = state.constants
    (; rho) = predictands

    if thermal_conductivity == 0.0
        return
    end

    phitheta .= 0.0

    mu_conduct = thermal_conductivity / uref / lref

    #-----------------------------------------
    #             Zonal fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        coef_t =
            mu_conduct *
            0.5 *
            (
                rhobar[i, j, k0] / rhobar[i, j, k] +
                rhobar[i + 1, j, k0] / rhobar[i + 1, j, k]
            )

        thetal = pbar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        thetar = pbar[i + 1, j, k] / (rho[i + 1, j, k] + rhobar[i + 1, j, k])

        thetad =
            0.5 * (
                pbar[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1]) +
                pbar[i + 1, j, k - 1] /
                (rho[i + 1, j, k - 1] + rhobar[i + 1, j, k - 1])
            )
        thetau =
            0.5 * (
                pbar[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1]) +
                pbar[i + 1, j, k + 1] /
                (rho[i + 1, j, k + 1] + rhobar[i + 1, j, k + 1])
            )

        dtht_dxi =
            0.5 * (jac[i, j, k] + jac[i + 1, j, k]) * (thetar - thetal) / dx +
            0.5 *
            (
                jac[i, j, k] * met[i, j, k, 1, 3] +
                jac[i + 1, j, k] * met[i + 1, j, k, 1, 3]
            ) *
            (thetau - thetad) / (2.0 * dz)

        phitheta[i, j, k, 1] = -coef_t * dtht_dxi
    end

    #-----------------------------------------
    #           Meridional fluxes
    #-----------------------------------------

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        coef_t =
            mu_conduct *
            0.5 *
            (
                rhobar[i, j, k0] / rhobar[i, j, k] +
                rhobar[i, j + 1, k0] / rhobar[i, j + 1, k]
            )

        thetab = pbar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        thetaf = pbar[i, j + 1, k] / (rho[i, j + 1, k] + rhobar[i, j + 1, k])

        thetad =
            0.5 * (
                pbar[i, j, k - 1] / (rho[i, j, k - 1] + rhobar[i, j, k - 1]) +
                pbar[i, j + 1, k - 1] /
                (rho[i, j + 1, k - 1] + rhobar[i, j + 1, k - 1])
            )
        thetau =
            0.5 * (
                pbar[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1]) +
                pbar[i, j + 1, k + 1] /
                (rho[i, j + 1, k + 1] + rhobar[i, j + 1, k + 1])
            )

        dtht_dyi =
            0.5 * (jac[i, j, k] + jac[i, j + 1, k]) * (thetaf - thetab) / dy +
            0.5 *
            (
                jac[i, j, k] * met[i, j, k, 2, 3] +
                jac[i, j + 1, k] * met[i, j + 1, k, 2, 3]
            ) *
            (thetau - thetad) / (2.0 * dz)

        phitheta[i, j, k, 2] = -coef_t * dtht_dyi
    end

    #-----------------------------------------
    #            Vertical fluxes
    #-----------------------------------------

    @ivy for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
        coef_t =
            mu_conduct * (
                jac[i, j, k + 1] * rhobar[i, j, 1] / rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, 1] / rhobar[i, j, k + 1]
            ) / (jac[i, j, k + 1] + jac[i, j, k])

        thetal =
            (
                jac[i - 1, j, k + 1] * pbar[i - 1, j, k] /
                (rho[i - 1, j, k] + rhobar[i - 1, j, k]) +
                jac[i - 1, j, k] * pbar[i - 1, j, k + 1] /
                (rho[i - 1, j, k + 1] + rhobar[i - 1, j, k + 1])
            ) / (jac[i - 1, j, k + 1] + jac[i - 1, j, k])

        thetar =
            (
                jac[i + 1, j, k + 1] * pbar[i + 1, j, k] /
                (rho[i + 1, j, k] + rhobar[i + 1, j, k]) +
                jac[i + 1, j, k] * pbar[i + 1, j, k + 1] /
                (rho[i + 1, j, k + 1] + rhobar[i + 1, j, k + 1])
            ) / (jac[i + 1, j, k + 1] + jac[i + 1, j, k])

        thetab =
            (
                jac[i, j - 1, k + 1] * pbar[i, j - 1, k] /
                (rho[i, j - 1, k] + rhobar[i, j - 1, k]) +
                jac[i, j - 1, k] * pbar[i, j - 1, k + 1] /
                (rho[i, j - 1, k + 1] + rhobar[i, j - 1, k + 1])
            ) / (jac[i, j - 1, k + 1] + jac[i, j - 1, k])

        thetaf =
            (
                jac[i, j + 1, k + 1] * pbar[i, j + 1, k] /
                (rho[i, j + 1, k] + rhobar[i, j + 1, k]) +
                jac[i, j + 1, k] * pbar[i, j + 1, k + 1] /
                (rho[i, j + 1, k + 1] + rhobar[i, j + 1, k + 1])
            ) / (jac[i, j + 1, k + 1] + jac[i, j + 1, k])

        thetad = pbar[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        thetau = pbar[i, j, k + 1] / (rho[i, j, k + 1] + rhobar[i, j, k + 1])

        dtht_dzi =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (
                (met[i, j, k, 1, 3] + met[i, j, k + 1, 1, 3]) *
                (thetar - thetal) / (2.0 * dx) +
                (met[i, j, k, 2, 3] + met[i, j, k + 1, 2, 3]) *
                (thetaf - thetab) / (2.0 * dy) +
                (met[i, j, k, 3, 3] + met[i, j, k + 1, 3, 3]) *
                (thetau - thetad) / dz
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        phitheta[i, j, k, 3] = -coef_t * dtht_dzi
    end

    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulence_scheme::NoTurbulence,
)
    return
end

function compute_fluxes!(
    state::State,
    predictands::Predictands,
    turbulence_scheme::AbstractTurbulence,
)

    # Get all necessary fields.
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; pbar) = state.atmosphere
    (; turbulencereconstructions, turbulencefluxes) = state.turbulence

    # Get old wind.
    (u0, v0, w0) = (predictands.u, predictands.v, predictands.w)

    @ivy for field in 1:fieldcount(TurbulencePredictands)
        chir = getfield(turbulencereconstructions, field)[2:end, :, :, 1, 1]
        chil = getfield(turbulencereconstructions, field)[:, :, :, 1, 2]
        fchi = getfield(turbulencefluxes, field)[:, :, :, 1]
        for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
            pedger =
                0.5 * (
                    jac[i, j, k] * pbar[i, j, k] +
                    jac[i + 1, j, k] * pbar[i + 1, j, k]
                )
            usurf = pedger * u0[i, j, k]

            fchi[i, j, k] = compute_flux(usurf, chil[i, j, k], chir[i, j, k])
        end

        chif = getfield(turbulencereconstructions, field)[:, 2:end, :, 2, 1]
        chib = getfield(turbulencereconstructions, field)[:, :, :, 2, 2]
        gchi = getfield(turbulencefluxes, field)[:, :, :, 2]
        for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
            pedgef =
                0.5 * (
                    jac[i, j, k] * pbar[i, j, k] +
                    jac[i, j + 1, k] * pbar[i, j + 1, k]
                )
            vsurf = pedgef * v0[i, j, k]

            gchi[i, j, k] = compute_flux(vsurf, chib[i, j, k], chif[i, j, k])
        end

        chiu = getfield(turbulencereconstructions, field)[:, :, 2:end, 3, 1]
        chid = getfield(turbulencereconstructions, field)[:, :, :, 3, 2]
        hchi = getfield(turbulencefluxes, field)[:, :, :, 3]
        for k in (k0 - 1):k1, j in j0:j1, i in i0:i1
            pedgeu =
                jac[i, j, k] *
                jac[i, j, k + 1] *
                (pbar[i, j, k] + pbar[i, j, k + 1]) /
                (jac[i, j, k] + jac[i, j, k + 1])
            wsurf = pedgeu * w0[i, j, k]

            hchi[i, j, k] = compute_flux(wsurf, chid[i, j, k], chiu[i, j, k])
        end
    end

    return
end