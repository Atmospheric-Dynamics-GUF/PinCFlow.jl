"""
```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::Rho)
```

Update the density if the atmosphere is not Boussinesq by dispatching to the appropriate method.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Boussinesq,
)
```

Return in Boussinesq mode (the density is constant).

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Union{PseudoIncompressible, Compressible},
)
```

Update the density with a Runge-Kutta step on the left-hand side of the equation (the right-hand side is zero).

The update is given by

```math
\\begin{align*}
    q^\\rho & \\rightarrow - \\frac{\\Delta t}{J} \\left(\\frac{\\mathcal{F}^{\\rho, \\widehat{x}}_{i + 1 / 2} - \\mathcal{F}^{\\rho, \\widehat{x}}_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho, \\widehat{y}}_{j + 1 / 2} - \\mathcal{F}^{\\rho, \\widehat{y}}_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}^{\\rho, \\widehat{z}}_{k + 1 / 2} - \\mathcal{F}^{\\rho, \\widehat{z}}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + \\alpha_\\mathrm{RK} q^\\rho,\\\\
    \\rho & \\rightarrow \\rho + \\beta_\\mathrm{RK} q^\\rho,
\\end{align*}
```

where ``\\Delta t`` is the time step given as input to this method.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::RhoP, side::LHS)
```

Update the density fluctuations with a Runge-Kutta step on the left-hand-side of the equation.

The update is given by

```math
\\begin{align*}
    q^{\\rho'} & \\rightarrow - \\frac{\\Delta t}{J} \\left(\\frac{\\mathcal{F}^{\\rho', \\widehat{x}}_{i + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{x}}_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho', \\widehat{y}}_{j + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{y}}_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}^{\\rho', \\widehat{z}}_{k + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{z}}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + \\alpha_\\mathrm{RK} q^{\\rho'},\\\\
    \\rho' & \\rightarrow \\rho' + \\beta_\\mathrm{RK} q^{\\rho'}
\\end{align*}
```

in Boussinesq/pseudo-incompressible mode and

```math
\\begin{align*}
    q^{\\rho'} & \\rightarrow \\Delta t \\left[- \\frac{1}{J} \\left(\\frac{\\mathcal{F}^{\\rho', \\widehat{x}}_{i + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{x}}_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho', \\widehat{y}}_{j + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{y}}_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}^{\\rho', \\widehat{z}}_{k + 1 / 2} - \\mathcal{F}^{\\rho', \\widehat{z}}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + \\frac{F^P}{\\overline{\\theta}}\\right] + \\alpha_\\mathrm{RK} q^{\\rho'},\\\\
    \\rho' & \\rightarrow \\rho' + \\beta_\\mathrm{RK} q^{\\rho'}
\\end{align*}
```

in compressible mode.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Explicit,
)
```

Update the density fluctuations with an explicit Euler step the on right-hand side of the equation, without the Rayleigh-damping term.

The update is given by

```math
\\rho' \\rightarrow - \\frac{\\rho}{g} \\left(b' - \\Delta t N^2 \\frac{\\overline{\\rho}}{\\rho} w\\right)
```

in Boussinesq/pseudo-incompressible mode and

```math
\\rho' \\rightarrow - \\frac{\\rho}{g} \\left[b' - \\Delta t N^2 \\frac{P / \\overline{\\theta}}{\\rho} \\left(\\frac{W_{k + 1 / 2}}{\\left(J P\\right)_{k + 1 / 2}}\\right)\\right]
```

in compressible mode, where ``b' = - g \\rho' / \\rho``.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
```

Update the density fluctuations with an implicit Euler step on the right-hand side of the equation.

The update is given by

```math
\\begin{align*}
    \\rho' & \\rightarrow - \\frac{\\rho}{g} \\left[1 + \\beta_\\mathrm{R} \\Delta t + \\frac{\\overline{\\rho}}{\\rho} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{- \\frac{\\overline{\\rho}}{\\rho} N^2 \\Delta t J \\left[\\widehat{w}_\\mathrm{old} + \\Delta t \\left(- \\left(c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}}\\right) + \\left(\\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)\\right)\\right] + \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right) b'\\right.\\\\
    & \\qquad \\quad + \\left.\\frac{\\overline{\\rho}}{\\rho} N^2 \\Delta t J \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right) \\left(G^{13} u + G^{23} v\\right)\\vphantom{\\left[\\left(\\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)\\right]}\\right\\},
\\end{align*}
```

in Boussinesq/pseudo-incompressible mode and

```math
\\begin{align*}
    \\rho' & \\rightarrow - \\frac{\\rho}{g} \\left[1 + \\beta_\\mathrm{R} \\Delta t + \\frac{P / \\overline{\\theta}}{\\rho} \\left(N \\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{- \\frac{P / \\overline{\\theta}}{\\rho} N^2 \\Delta t J \\left[\\left(\\frac{\\widehat{W}_{\\mathrm{old}, k + 1 / 2}}{\\left(J P\\right)_{k + 1 / 2}}\\right) + \\Delta t \\left(- \\left(c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}}\\right) + \\left(\\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)\\right)\\right]\\right.\\\\
    & \\qquad \\quad + \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right) b' + \\frac{P / \\overline{\\theta}}{\\rho} N^2 \\Delta t J \\left(1 + \\beta_\\mathrm{R} \\Delta t\\right)\\\\
    & \\qquad \\quad \\times \\left.\\left[G^{13} \\left(\\frac{U_{i + 1 / 2}}{\\left(J P\\right)_{i + 1 / 2}}\\right) + G^{23} \\left(\\frac{V_{j + 1 / 2}}{\\left(J P\\right)_{j + 1 / 2}}\\right)\\right]\\right\\},
\\end{align*}
```

in compressible mode, where ``\\widehat{w}_\\mathrm{old}`` is the transformed vertical wind stored in `state.variables.backups`.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::U, side::LHS)
```

Update the zonal momentum with a Runge-Kutta step on the left-hand side of the equation.

The update is given by

```math
\\begin{align*}
    q^{\\rho u}_{i + 1 / 2} & \\rightarrow \\Delta t \\left[- \\frac{1}{J_{i + 1 / 2}} \\left(\\frac{\\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} - \\mathcal{F}^{\\rho u, \\widehat{x}}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} - \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j - 1 / 2}}{\\Delta \\widehat{y}}\\right.\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad + \\left.\\left.\\frac{\\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + f \\left(\\rho_\\mathrm{old} v\\right)_{i + 1 / 2}\\right] + \\alpha_\\mathrm{RK} q^{\\rho u}_{i + 1 / 2},\\\\
    u_{i + 1 / 2} & \\rightarrow \\rho_{i + 1 / 2}^{- 1} \\left(\\rho_{\\mathrm{old}, i + 1 / 2} u_{i + 1 / 2} + \\beta_\\mathrm{RK} q^{\\rho u}_{i + 1 / 2}\\right),
\\end{align*}
```

where ``\\rho_\\mathrm{old}`` is the density stored in `state.variables.backups`.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Explicit,
)
```

Update the zonal wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

The update is given by

```math
u_{i + 1 / 2} \\rightarrow u_{i + 1 / 2} + \\Delta t \\left(- c_p \\frac{P_{i + 1 / 2}}{\\rho_{i + 1 / 2}} \\mathcal{P}_{i + 1 / 2}^{\\rho u} + \\frac{F_{i + 1 / 2}^{\\rho u}}{\\rho_{i + 1 / 2}}\\right)
```

in Boussinesq/pseudo-incompressible mode and

```math
U_{i + 1 / 2} \\rightarrow U_{i + 1 / 2} + \\Delta t \\left(J P\\right)_{i + 1 / 2} \\left(- c_p \\frac{P_{i + 1 / 2}}{\\rho_{i + 1 / 2}} \\mathcal{P}_{i + 1 / 2}^{\\rho u} + \\frac{F_{i + 1 / 2}^{\\rho u}}{\\rho_{i + 1 / 2}}\\right)
```

in compressible mode.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
```

Update the zonal wind with an implicit Euler step on the right-hand side of the equation.

The update is given by

```math
u_{i + 1 / 2} \\rightarrow \\left(1 + \\beta_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\left[u_{i + 1 / 2} + \\Delta t \\left(- c_p \\frac{P_{i + 1 / 2}}{\\rho_{i + 1 / 2}} \\mathcal{P}_{i + 1 / 2}^{\\rho u} + \\frac{F_{i + 1 / 2}^{\\rho u}}{\\rho_{i + 1 / 2}}\\right)\\right]
```

in Boussinesq/pseudo-incompressible mode and

```math
U_{i + 1 / 2} \\rightarrow \\left(1 + \\beta_{\\mathrm{R}, i + 1 / 2} \\Delta t\\right)^{- 1} \\left[U_{i + 1 / 2} + \\Delta t \\left(J P\\right)_{i + 1 / 2} \\left(- c_p \\frac{P_{i + 1 / 2}}{\\rho_{i + 1 / 2}} \\mathcal{P}_{i + 1 / 2}^{\\rho u} + \\frac{F_{i + 1 / 2}^{\\rho u}}{\\rho_{i + 1 / 2}}\\right)\\right]
```

in compressible mode.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::V, side::LHS)
```

Update the meridional momentum with a Runge-Kutta step on the left-hand side of the equation.

The update is given by

```math
\\begin{align*}
    q^{\\rho v}_{j + 1 / 2} & \\rightarrow \\Delta t \\left[- \\frac{1}{J_{j + 1 / 2}} \\left(\\frac{\\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} - \\mathcal{F}^{\\rho v, \\widehat{x}}_{i - 1 / 2, j + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} - \\mathcal{F}^{\\rho v, \\widehat{y}}}{\\Delta \\widehat{y}}\\right.\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad + \\left.\\left.\\frac{\\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k - 1 / 2}}{\\Delta \\widehat{z}}\\right) - f \\left(\\rho_\\mathrm{old} u_\\mathrm{old}\\right)_{j + 1 / 2}\\right] + \\alpha_\\mathrm{RK} q^{\\rho v}_{j + 1 / 2},\\\\
    v_{j + 1 / 2} & \\rightarrow \\rho_{j + 1 / 2}^{- 1} \\left(\\rho_{\\mathrm{old}, j + 1 / 2} v_{j + 1 / 2} + \\beta_\\mathrm{RK} q^{\\rho v}_{j + 1 / 2}\\right),
\\end{align*}
```

where ``\\rho_\\mathrm{old}`` and ``u_{\\mathrm{old}, i + 1 / 2}`` are the density and zonal wind stored in `state.variables.backups`.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Explicit,
)
```

Update the meridional wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

The update is given by

```math
v_{i + 1 / 2} \\rightarrow v_{j + 1 / 2} + \\Delta t \\left(- c_p \\frac{P_{j + 1 / 2}}{\\rho_{j + 1 / 2}} \\mathcal{P}_{j + 1 / 2}^{\\rho v} + \\frac{F_{j + 1 / 2}^{\\rho v}}{\\rho_{j + 1 / 2}}\\right)
```

in Boussinesq/pseudo-incompressible mode and

```math
V_{j + 1 / 2} \\rightarrow V_{j + 1 / 2} + \\Delta t \\left(J P\\right)_{j + 1 / 2} \\left(- c_p \\frac{P_{j + 1 / 2}}{\\rho_{j + 1 / 2}} \\mathcal{P}_{j + 1 / 2}^{\\rho v} + \\frac{F_{j + 1 / 2}^{\\rho v}}{\\rho_{j + 1 / 2}}\\right)
```

in compressible mode.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
```

Update the meridional wind with an implicit Euler step on the right-hand side of the equation.

The update is given by

```math
v_{j + 1 / 2} \\rightarrow \\left(1 + \\beta_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\left[v_{j + 1 / 2} + \\Delta t \\left(- c_p \\frac{P_{j + 1 / 2}}{\\rho_{j + 1 / 2}} \\mathcal{P}_{j + 1 / 2}^{\\rho v} + \\frac{F_{j + 1 / 2}^{\\rho v}}{\\rho_{j + 1 / 2}}\\right)\\right]
```

in Boussinesq/pseudo-incompressible mode and

```math
V_{j + 1 / 2} \\rightarrow \\left(1 + \\beta_{\\mathrm{R}, j + 1 / 2} \\Delta t\\right)^{- 1} \\left[V_{j + 1 / 2} + \\Delta t \\left(J P\\right)_{j + 1 / 2} \\left(- c_p \\frac{P_{j + 1 / 2}}{\\rho_{j + 1 / 2}} \\mathcal{P}_{j + 1 / 2}^{\\rho v} + \\frac{F_{j + 1 / 2}^{\\rho v}}{\\rho_{j + 1 / 2}}\\right)\\right]
```

in compressible mode.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::W, side::LHS)
```

Update the transformed vertical momentum with a Runge-Kutta step on the left-hand side of the equation.

The update is given by

```math
\\begin{align*}
    q^{\\rho \\widehat{w}}_{k + 1 / 2} & \\rightarrow \\Delta t \\left\\{- \\left[G^{13} \\left(\\frac{1}{J_{i + 1 / 2}} \\left(\\frac{\\mathcal{F}^{\\rho u, \\widehat{x}}_{i + 1} - \\mathcal{F}^{\\rho u, \\widehat{x}}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j + 1 / 2} - \\mathcal{F}^{\\rho u, \\widehat{y}}_{i + 1 / 2, j - 1 / 2}}{\\Delta \\widehat{y}}\\right.\\right.\\right.\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad \\qquad \\qquad + \\left.\\left.\\left.\\frac{\\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho u, \\widehat{z}}_{i + 1 / 2, k - 1 / 2}}{\\Delta \\widehat{z}}\\right)\\right)\\right]_{k + 1 / 2}\\\\
    & \\qquad \\qquad - \\left[G^{23} \\left(\\frac{1}{J_{j + 1 / 2}} \\left(\\frac{\\mathcal{F}^{\\rho v, \\widehat{x}}_{i + 1 / 2, j + 1 / 2} - \\mathcal{F}^{\\rho v, \\widehat{x}}_{i - 1 / 2, j + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho v, \\widehat{y}}_{j + 1} - \\mathcal{F}^{\\rho v, \\widehat{y}}}{\\Delta \\widehat{y}}\\right.\\right.\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad \\qquad \\qquad + \\left.\\left.\\left.\\frac{\\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho v, \\widehat{z}}_{j + 1 / 2, k - 1 / 2}}{\\Delta \\widehat{z}}\\right)\\right)\\right]_{k + 1 / 2}\\\\
    & \\qquad \\qquad - \\frac{1}{J_{k + 1 / 2}^2} \\left(\\frac{\\mathcal{F}^{\\rho w, \\widehat{x}}_{i + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho w, \\widehat{x}}_{i - 1 / 2, k + 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho w, \\widehat{y}}_{j + 1 / 2, k + 1 / 2} - \\mathcal{F}^{\\rho w, \\widehat{y}}_{j - 1 / 2, k + 1 / 2}}{\\Delta \\widehat{y}}\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad \\quad + \\left.\\frac{\\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} - \\mathcal{F}^{\\rho w, \\widehat{z}}}{\\Delta \\widehat{z}}\\right)\\\\
    & \\qquad \\qquad + \\left.G^{13} f \\left(\\rho_\\mathrm{old} v_\\mathrm{old}\\right)_{k + 1 / 2} - G^{23} f \\left(\\rho_\\mathrm{old} u_\\mathrm{old}\\right)_{k + 1 / 2}\\vphantom{- \\frac{1}{J^2} \\left(\\frac{\\mathcal{F}^{\\rho w, \\widehat{z}}_{k + 1} - \\mathcal{F}^{\\rho w, \\widehat{z}}}{\\Delta \\widehat{z}}\\right)}\\right\\} + \\alpha_\\mathrm{RK} q^{\\rho \\widehat{w}}_{k + 1 / 2},\\\\
    \\widehat{w}_{k + 1 / 2} & \\rightarrow \\rho_{k + 1 / 2}^{- 1} \\left(\\rho_{\\mathrm{old}, k + 1 / 2} \\widehat{w}_{k + 1 / 2} + \\beta_\\mathrm{RK} q^{\\rho \\widehat{w}}_{k + 1 / 2}\\right),
\\end{align*}
```

where ``\\rho_\\mathrm{old}``, ``u_{\\mathrm{old}, i + 1 / 2}`` and ``v_{\\mathrm{old}, j + 1 / 2}`` are the density, zonal wind and meridional wind stored in `state.variables.backups`.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Explicit,
)
```

Update the transformed vertical wind with an explicit Euler step on the right-hand side of the equation, without the Rayleigh-damping term.

The update is given by

```math
\\widehat{w}_{k + 1 / 2} \\rightarrow \\widehat{w}_{k + 1 / 2} + \\Delta t \\left[- c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}} + \\left(\\frac{b'_\\mathrm{old}}{J}\\right)_{k + 1 / 2} + \\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right]
```

in Boussinesq/pseudo-incompressible mode and

```math
\\widehat{W}_{k + 1 / 2} \\rightarrow \\widehat{W}_{k + 1 / 2} + \\Delta t \\left(J P\\right)_{k + 1 / 2} \\left[- c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}} + \\left(\\frac{b'_\\mathrm{old}}{J}\\right)_{k + 1 / 2} + \\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right]
```

in compressible mode, where ``b'_\\mathrm{old} = - g \\rho'_\\mathrm{old} / \\rho``, with ``\\rho'_\\mathrm{old}`` being the density fluctuations stored in `state.variables.backups`.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
```

Update the transformed vertical wind with an implicit Euler step on the right-hand side of the equation.

The update is given by

```math
\\begin{align*}
    \\widehat{w}_{k + 1 / 2} & \\rightarrow \\left[1 + \\beta_{\\mathrm{R}, k + 1 / 2} \\Delta t + \\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} N^2_{k + 1 / 2} \\left(\\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{\\widehat{w}_{k + 1 / 2} + \\Delta t \\left(- c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}} + \\left(\\frac{b'}{J}\\right)_{k + 1 / 2} + \\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)\\right.\\\\
    & \\qquad \\quad + \\left.\\frac{\\overline{\\rho}_{k + 1 / 2}}{\\rho_{k + 1 / 2}} N^2_{k + 1 / 2} \\left(\\Delta t\\right)^2 \\left[\\left(G^{13} u\\right)_{k + 1 / 2} + \\left(G^{2 3} v\\right)_{k + 1 / 2}\\right]\\vphantom{\\left(\\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)}\\right\\}
\\end{align*}
```

in Boussinesq/pseudo-incompressible mode and

```math
\\begin{align*}
    \\widehat{W}_{k + 1 / 2} & \\rightarrow \\left[1 + \\beta_{\\mathrm{R}, k + 1 / 2} \\Delta t + \\frac{\\left(P / \\overline{\\theta}\\right)_{k + 1 / 2}}{\\rho_{k + 1 / 2}} N^2_{k + 1 / 2} \\left(\\Delta t\\right)^2\\right]^{- 1}\\\\
    & \\quad \\times \\left\\{\\widehat{W}_{k + 1 / 2} + \\Delta t \\left(J P\\right)_{k + 1 / 2} \\left(- c_p \\frac{P_{k + 1 / 2}}{\\rho_{k + 1 / 2}} \\mathcal{P}_{k + 1 / 2}^{\\rho \\widehat{w}} + \\left(\\frac{b'}{J}\\right)_{k + 1 / 2} + \\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)\\right.\\\\
    & \\qquad \\quad + \\left(J P\\right)_{k + 1 / 2} \\frac{\\left(P / \\overline{\\theta}\\right)_{k + 1 / 2}}{\\rho_{k + 1 / 2}} N^2_{k + 1 / 2} \\left(\\Delta t\\right)^2\\\\
    & \\qquad \\quad \\times \\left.\\left[\\left(G^{13} \\left(\\frac{U_{i + 1 / 2}}{\\left(J P\\right)_{i + 1 / 2}}\\right)\\right)_{k + 1 / 2} + \\left(G^{2 3} \\left(\\frac{V_{j + 1 / 2}}{\\left(J P\\right)_{j + 1 / 2}}\\right)\\right)_{k + 1 / 2}\\right]\\vphantom{\\left(\\frac{F_{k + 1 / 2}^{\\rho \\widehat{w}}}{\\rho_{k + 1 / 2}}\\right)}\\right\\}
\\end{align*}
```

in compressible mode.

```julia
update!(state::State, dt::AbstractFloat, variable::PiP)
```

Update the Exner-pressure if the atmosphere is compressible by dispatching to the appropriate method.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    variable::PiP,
    model::Union{Boussinesq, PseudoIncompressible},
)
```

Return in non-compressible modes.

```julia
update!(state::State, dt::AbstractFloat, variable::PiP, model::Compressible)
```

Update the Exner-pressure such that it is synchronized with the updated mass-weighted potential temperature.

The update is given by

```math
\\begin{align*}
    \\pi' & \\rightarrow \\pi' + \\Delta t \\left(\\frac{\\partial \\pi'}{\\partial P}\\right) \\left[- \\frac{1}{J} \\left(\\frac{U_{\\mathrm{old}, i + 1 / 2} - U_{\\mathrm{old}, i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{V_{\\mathrm{old}, j + 1 / 2} - V_{\\mathrm{old}, j - 1 / 2}}{\\Delta \\widehat{y}}\\right.\\right.\\\\
    & \\qquad \\qquad \\qquad \\qquad \\qquad \\qquad + \\left.\\left.\\frac{\\widehat{W}_{\\mathrm{old}, k + 1 / 2} - \\widehat{W}_{\\mathrm{old}, k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + F^P\\right],
\\end{align*}
```

where ``U_{\\mathrm{old}, i + 1 / 2}``, ``V_{\\mathrm{old}, j + 1 / 2}`` and ``\\widehat{W}_{\\mathrm{old}, k + 1 / 2}`` are the transformed wind components (including the factor ``J P``) stored in `state.variables.backups`.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, variable::P)
```

Update the mass-weighted potential temperature if the atmosphere is compressible by dispatching to the appropriate method.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Union{Boussinesq, PseudoIncompressible},
)
```

Return in non-compressible modes.

```julia
update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Compressible,
)
```

Update the mass-weighted potential temperature with a Runge-Kutta step on the left-hand side of the equation (the right-hand side is zero).

The update is given by

```math
\\begin{align*}
    q^P & \\rightarrow \\Delta t \\left[- \\frac{1}{J} \\left(\\frac{\\mathcal{F}^{P, \\widehat{x}}_{i + 1 / 2} - \\mathcal{F}^{P, \\widehat{x}}_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{P, \\widehat{y}}_{j + 1 / 2} - \\mathcal{F}^{P, \\widehat{y}}_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}^{P, \\widehat{z}}_{k + 1 / 2} - \\mathcal{F}^{P, \\widehat{z}}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + F^P\\right] + \\alpha_\\mathrm{RK} q^P,\\\\
    P & \\rightarrow P + \\beta_\\mathrm{RK} q^P.
\\end{align*}
```

```julia
update!(state::State, dt::AbstractFloat, m::Integer, tracer_setup::NoTracer)
```

Return for configurations without tracer transport.

```julia
update!(state::State, dt::AbstractFloat, m::Integer, tracer_setup::TracerOn)
```

Update the tracers with a Runge-Kutta step on the left-hand sides of the equations with WKB right-hand side terms according to namelists configuration.

The update is given by

```math
\\begin{align*}
    q^{\\rho \\chi} & \\rightarrow \\Delta t \\left[- \\frac{1}{J} \\left(\\frac{\\mathcal{F}^{\\rho \\chi, \\widehat{x}}_{i + 1 / 2} - \\mathcal{F}^{\\rho \\chi, \\widehat{x}}_{i - 1 / 2}}{\\Delta \\widehat{x}} + \\frac{\\mathcal{F}^{\\rho \\chi, \\widehat{y}}_{j + 1 / 2} - \\mathcal{F}^{\\rho \\chi, \\widehat{y}}_{j - 1 / 2}}{\\Delta \\widehat{y}} + \\frac{\\mathcal{F}^{\\rho \\chi, \\widehat{z}}_{k + 1 / 2} - \\mathcal{F}^{\\rho \\chi, \\widehat{z}}_{k - 1 / 2}}{\\Delta \\widehat{z}}\\right) + F^{\\rho \\chi}\\right] + \\alpha_\\mathrm{RK} q^{\\rho \\chi},\\\\
    \\left(\\rho \\chi\\right) & \\rightarrow \\left(\\rho \\chi\\right) + \\beta_\\mathrm{RK} q^{\\rho \\chi}.
\\end{align*}
```

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `m`: Runge-Kutta-stage index.

  - `variable`: Variable to update.

  - `model`: Dynamic equations.

  - `side`: Side of the equation.

  - `integration`: Type of the Euler step.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

  - `tracer_setup`: General tracer-transport configuration.

# See also

  - [`PinCFlow.Update.compute_volume_force`](@ref)

  - [`PinCFlow.Update.compute_compressible_wind_factor`](@ref)

  - [`PinCFlow.Update.compute_vertical_wind`](@ref)

  - [`PinCFlow.Update.compute_buoyancy_factor`](@ref)

  - [`PinCFlow.Update.compute_pressure_gradient`](@ref)

  - [`PinCFlow.Update.transform`](@ref)

  - [`PinCFlow.Update.conductive_heating`](@ref)
"""
function update! end

function update!(state::State, dt::AbstractFloat, m::Integer, variable::Rho)
    (; model) = state.namelists.setting
    update!(state, dt, m, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Boussinesq,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::Rho,
    model::Union{PseudoIncompressible, Compressible},
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; drho) = state.variables.increments
    (; phirho) = state.variables.fluxes
    (; rho) = state.variables.predictands

    if m == 1
        drho .= 0.0
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirho[i - 1, j, k, 1]
        fr = phirho[i, j, k, 1]
        gb = phirho[i, j - 1, k, 2]
        gf = phirho[i, j, k, 2]
        hd = phirho[i, j, k - 1, 3]
        hu = phirho[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        f = -fluxdiff

        drho[i, j, k] = dt * f + alphark[m] * drho[i, j, k]
        rho[i, j, k] += betark[m] * drho[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::RhoP,
    side::LHS,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; thetabar) = state.atmosphere
    (; alphark, betark) = state.time
    (; drhop) = state.variables.increments
    (; phirhop) = state.variables.fluxes
    (; rhop) = state.variables.predictands

    if m == 1
        drhop .= 0.0
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phirhop[i - 1, j, k, 1]
        fr = phirhop[i, j, k, 1]
        gb = phirhop[i, j - 1, k, 2]
        gf = phirhop[i, j, k, 2]
        hd = phirhop[i, j, k - 1, 3]
        hu = phirhop[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, i, j, k, P())

        f = -fluxdiff - heating / thetabar[i, j, k]

        drhop[i, j, k] = dt * f + alphark[m] * drhop[i, j, k]
        rhop[i, j, k] += betark[m] * drhop[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Explicit,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; g_ndim) = state.constants
    (; rhobar, n2) = state.atmosphere
    (; predictands) = state.variables
    (; rho, rhop) = predictands

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        jpu = compute_compressible_wind_factor(state, i, j, k, W())
        jpd = compute_compressible_wind_factor(state, i, j, k - 1, W())
        wvrt =
            0.5 * (
                compute_vertical_wind(i, j, k, state) / jpu +
                compute_vertical_wind(i, j, k - 1, state) / jpd
            )

        buoy = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        fb = compute_buoyancy_factor(state, i, j, k, RhoP())
        buoy -= dt * fb * n2[i, j, k] * wvrt

        rhop[i, j, k] = -buoy * (rho[i, j, k] + rhobar[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::RhoP,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; zz_size, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; betar) = state.sponge
    (; g_ndim) = state.constants
    (; rhobar, n2) = state.atmosphere
    (; rho, rhop, u, v, pip) = state.variables.predictands
    (; wold) = state.variables.backups

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        rhoc = rho[i, j, k] + rhobar[i, j, k]
        rhoedgeu =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhoedgeu +=
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhoedged =
            (
                jac[i, j, k - 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        rhoedged +=
            (
                jac[i, j, k - 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])

        jpedgeu = compute_compressible_wind_factor(state, i, j, k, W())
        jpedged = compute_compressible_wind_factor(state, i, j, k - 1, W())
        w = 0.5 * (wold[i, j, k] / jpedgeu + wold[i, j, k - 1] / jpedged)

        lower_gradient = compute_pressure_gradient(state, pip, i, j, k - 1, W())
        lower_force = compute_volume_force(state, i, j, k - 1, W())
        upper_gradient = compute_pressure_gradient(state, pip, i, j, k, W())
        upper_force = compute_volume_force(state, i, j, k, W())

        if ko + k == k0
            lower_gradient = 0.0
            lower_force = 0.0
        elseif ko + k == zz_size - nbz
            upper_gradient = 0.0
            upper_force = 0.0
        end

        gradient = 0.5 * (lower_gradient + upper_gradient)
        force = 0.5 * (lower_force / rhoedged + upper_force / rhoedgeu) * rhoc

        factor = 1.0

        factor += dt * betar[i, j, k] * rayleigh_factor

        b = -g_ndim * rhop[i, j, k] / (rho[i, j, k] + rhobar[i, j, k])
        jpedger = compute_compressible_wind_factor(state, i, j, k, U())
        jpedgel = compute_compressible_wind_factor(state, i - 1, j, k, U())
        jpedgef = compute_compressible_wind_factor(state, i, j, k, V())
        jpedgeb = compute_compressible_wind_factor(state, i, j - 1, k, V())
        fb = compute_buoyancy_factor(state, i, j, k, RhoP())
        b =
            1.0 / (factor + fb * n2[i, j, k] * dt^2.0) * (
                -fb *
                n2[i, j, k] *
                dt *
                jac[i, j, k] *
                (w + dt * (-gradient + force / rhoc)) +
                factor * b +
                fb *
                n2[i, j, k] *
                dt *
                jac[i, j, k] *
                factor *
                0.5 *
                (
                    met[i, j, k, 1, 3] *
                    (u[i, j, k] / jpedger + u[i - 1, j, k] / jpedgel) +
                    met[i, j, k, 2, 3] *
                    (v[i, j, k] / jpedgef + v[i, j - 1, k] / jpedgeb)
                )
            )

        rhop[i, j, k] = -b * (rho[i, j, k] + rhobar[i, j, k]) / g_ndim
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::U,
    side::LHS,
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; alphark, betark) = state.time
    (; tref) = state.constants
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhobar) = state.atmosphere
    (; du) = state.variables.increments
    (; phiu) = state.variables.fluxes
    (; rhoold, uold) = state.variables.backups
    (; rho, u, v) = state.variables.predictands

    fc = coriolis_frequency * tref

    if m == 1
        du .= 0.0
    end

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1

        # Compute zonal momentum flux divergence.
        fr = phiu[i, j, k, 1]
        fl = phiu[i - 1, j, k, 1]
        gf = phiu[i, j, k, 2]
        gb = phiu[i, j - 1, k, 2]
        hu = phiu[i, j, k, 3]
        hd = phiu[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust zonal momentum flux divergence.
        jacedger = 0.5 * (jac[i, j, k] + jac[i + 1, j, k])
        fluxdiff /= jacedger

        # Explicit integration of Coriolis force in TFC.
        uold[i, j, k] = u[i, j, k]
        if k == k1 && ko + nzz != zz_size
            uold[i, j, k + 1] = u[i, j, k + 1]
        end
        vc = 0.5 * (v[i, j, k] + v[i, j - 1, k])
        vr = 0.5 * (v[i + 1, j, k] + v[i + 1, j - 1, k])
        volforce =
            0.5 *
            fc *
            (
                (rhoold[i, j, k] + rhobar[i, j, k]) * vc +
                (rhoold[i + 1, j, k] + rhobar[i + 1, j, k]) * vr
            )

        # Compute force.
        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i + 1, j, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhobaredger = 0.5 * (rhobar[i, j, k] + rhobar[i + 1, j, k])
        rhom_1 += rhobaredger
        rhom += rhobaredger

        # Set velocity and momentum at previous time.
        um_1 = u[i, j, k]
        momm_1 = rhom_1 * um_1

        # Compute tendency.
        du[i, j, k] = dt * force + alphark[m] * du[i, j, k]

        # Update momentum.
        momm = momm_1 + betark[m] * du[i, j, k]

        # Update wind.
        uast = momm / rhom
        u[i, j, k] = uast
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Explicit,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; rho, u, pip) = state.variables.predictands

    @ivy for k in k0:k1, j in j0:j1, i in (i0 - 1):i1
        rhoedger = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhobaredger = 0.5 * (rhobar[i, j, k] + rhobar[i + 1, j, k])
        rhoedger += rhobaredger

        gradient = compute_pressure_gradient(state, pip, i, j, k, U())

        force = compute_volume_force(state, i, j, k, U())

        jpedger = compute_compressible_wind_factor(state, i, j, k, U())

        u[i, j, k] += dt * (-gradient + force / rhoedger) * jpedger
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::U,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
    (; damp_horizontal_wind_on_rhs) = state.namelists.sponge
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; betar) = state.sponge
    (; rho, u, pip) = state.variables.predictands

    kmin = k0
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in j0:j1, i in (i0 - 1):i1
        rhoedger = 0.5 * (rho[i, j, k] + rho[i + 1, j, k])
        rhobaredger = 0.5 * (rhobar[i, j, k] + rhobar[i + 1, j, k])
        rhoedger += rhobaredger

        gradient = compute_pressure_gradient(state, pip, i, j, k, U())

        force = compute_volume_force(state, i, j, k, U())

        factor = 1.0

        if damp_horizontal_wind_on_rhs
            factor +=
                dt *
                0.5 *
                (betar[i, j, k] + betar[i + 1, j, k]) *
                rayleigh_factor
        end

        jpedger = compute_compressible_wind_factor(state, i, j, k, U())

        u[i, j, k] =
            1.0 / factor *
            (u[i, j, k] + dt * (-gradient + force / rhoedger) * jpedger)
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::V,
    side::LHS,
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; alphark, betark) = state.time
    (; tref) = state.constants
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; rhobar) = state.atmosphere
    (; dv) = state.variables.increments
    (; phiv) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, v) = state.variables.predictands

    fc = coriolis_frequency * tref

    if m == 1
        dv .= 0.0
    end

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1

        # Compute meridional momentum flux divergence.
        fr = phiv[i, j, k, 1]
        fl = phiv[i - 1, j, k, 1]
        gf = phiv[i, j, k, 2]
        gb = phiv[i, j - 1, k, 2]
        hu = phiv[i, j, k, 3]
        hd = phiv[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust meridional momentum flux divergence.
        jacedgef = 0.5 * (jac[i, j, k] + jac[i, j + 1, k])
        fluxdiff /= jacedgef

        # Explicit integration of Coriolis force in TFC.
        vold[i, j, k] = v[i, j, k]
        if k == k1 && ko + nzz != zz_size
            vold[i, j, k + 1] = v[i, j, k + 1]
        end
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uf = 0.5 * (uold[i, j + 1, k] + uold[i - 1, j + 1, k])

        volforce =
            -0.5 *
            fc *
            (
                (rhoold[i, j, k] + rhobar[i, j, k]) * uc +
                (rhoold[i, j + 1, k] + rhobar[i, j + 1, k]) * uf
            )

        force = -fluxdiff + volforce

        # Interpolate density.
        rhom_1 = 0.5 * (rhoold[i, j, k] + rhoold[i, j + 1, k])
        rhom = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhobaredgef = 0.5 * (rhobar[i, j, k] + rhobar[i, j + 1, k])
        rhom_1 += rhobaredgef
        rhom += rhobaredgef

        vm_1 = v[i, j, k]
        momm_1 = rhom_1 * vm_1

        dv[i, j, k] = dt * force + alphark[m] * dv[i, j, k]

        momm = momm_1 + betark[m] * dv[i, j, k]

        # Update wind.
        vast = momm / rhom
        v[i, j, k] = vast
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Explicit,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; rho, v, pip) = state.variables.predictands

    @ivy for k in k0:k1, j in (j0 - 1):j1, i in i0:i1
        rhoedgef = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhobaredgef = 0.5 * (rhobar[i, j, k] + rhobar[i, j + 1, k])
        rhoedgef += rhobaredgef

        gradient = compute_pressure_gradient(state, pip, i, j, k, V())

        force = compute_volume_force(state, i, j, k, V())

        jpedgef = compute_compressible_wind_factor(state, i, j, k, V())

        v[i, j, k] += dt * (-gradient + force / rhoedgef) * jpedgef
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::V,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
    (; damp_horizontal_wind_on_rhs) = state.namelists.sponge
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; rhobar) = state.atmosphere
    (; betar) = state.sponge
    (; rho, v, pip) = state.variables.predictands

    kmin = k0
    kmax = ko + nzz == zz_size ? k1 : k1 + 1

    @ivy for k in kmin:kmax, j in (j0 - 1):j1, i in i0:i1
        rhoedgef = 0.5 * (rho[i, j, k] + rho[i, j + 1, k])
        rhobaredgef = 0.5 * (rhobar[i, j, k] + rhobar[i, j + 1, k])
        rhoedgef += rhobaredgef

        gradient = compute_pressure_gradient(state, pip, i, j, k, V())

        force = compute_volume_force(state, i, j, k, V())

        factor = 1.0

        if damp_horizontal_wind_on_rhs
            factor +=
                dt *
                0.5 *
                (betar[i, j, k] + betar[i, j + 1, k]) *
                rayleigh_factor
        end

        jpedgef = compute_compressible_wind_factor(state, i, j, k, V())

        v[i, j, k] =
            1.0 / factor *
            (v[i, j, k] + dt * (-gradient + force / rhoedgef) * jpedgef)
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::W,
    side::LHS,
)
    (; coriolis_frequency) = state.namelists.atmosphere
    (; alphark, betark) = state.time
    (; tref) = state.constants
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; grid) = state
    (; dx, dy, dz, jac, met) = grid
    (; rhobar) = state.atmosphere
    (; dw) = state.variables.increments
    (; phiu, phiv, phiw) = state.variables.fluxes
    (; rhoold, uold, vold) = state.variables.backups
    (; rho, w) = state.variables.predictands

    fc = coriolis_frequency * tref

    # Initialize fields for transformation of momentum flux divergence.
    (fluxdiffu, fluxdiffv) = (zeros(2, 2) for i in 1:2)

    if m == 1
        dw .= 0.0
    end

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 - 1 : k1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        # Compute vertical momentum flux divergence.
        fr = phiw[i, j, k, 1]
        fl = phiw[i - 1, j, k, 1]
        gf = phiw[i, j, k, 2]
        gb = phiw[i, j - 1, k, 2]
        hu = phiw[i, j, k, 3]
        hd = phiw[i, j, k - 1, 3]
        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz

        # Adjust Cartesian vertical momentum flux divergence.
        jacedgeu =
            2.0 * jac[i, j, k] * jac[i, j, k + 1] /
            (jac[i, j, k] + jac[i, j, k + 1])
        fluxdiff /= jacedgeu

        # Compute zonal momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiu[i - ll, j, k + mm, 1]
            fl = phiu[i - 1 - ll, j, k + mm, 1]
            gf = phiu[i - ll, j, k + mm, 2]
            gb = phiu[i - ll, j - 1, k + mm, 2]
            hu = phiu[i - ll, j, k + mm, 3]
            hd = phiu[i - ll, j, k - 1 + mm, 3]
            fluxdiffu[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedger =
                0.5 * (jac[i - ll, j, k + mm] + jac[i + 1 - ll, j, k + mm])
            fluxdiffu[ll + 1, mm + 1] /= jacedger
        end

        # Compute meridional momentum flux divergences.
        for ll in 0:1, mm in 0:1
            fr = phiv[i, j - ll, k + mm, 1]
            fl = phiv[i - 1, j - ll, k + mm, 1]
            gf = phiv[i, j - ll, k + mm, 2]
            gb = phiv[i, j - 1 - ll, k + mm, 2]
            hu = phiv[i, j - ll, k + mm, 3]
            hd = phiv[i, j - ll, k - 1 + mm, 3]
            fluxdiffv[ll + 1, mm + 1] =
                (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            jacedgef =
                0.5 * (jac[i, j - ll, k + mm] + jac[i, j + 1 - ll, k + mm])
            fluxdiffv[ll + 1, mm + 1] /= jacedgef
        end

        # Compute transformed vertical momentum flux divergence.
        fluxdiff = transform(
            i,
            j,
            k,
            fluxdiffu[1, 1],
            fluxdiffu[1, 2],
            fluxdiffu[2, 1],
            fluxdiffu[2, 2],
            fluxdiffv[1, 1],
            fluxdiffv[1, 2],
            fluxdiffv[2, 1],
            fluxdiffv[2, 2],
            fluxdiff,
            Transformed(),
            state,
        )

        # Explicit integration of Coriolis force in TFC.
        vc = 0.5 * (vold[i, j, k] + vold[i, j - 1, k])
        vu = 0.5 * (vold[i, j, k + 1] + vold[i, j - 1, k + 1])
        uc = 0.5 * (uold[i, j, k] + uold[i - 1, j, k])
        uu = 0.5 * (uold[i, j, k + 1] + uold[i - 1, j, k + 1])

        volforce =
            fc * (
                jac[i, j, k + 1] *
                met[i, j, k, 1, 3] *
                (rhoold[i, j, k] + rhobar[i, j, k]) *
                vc +
                jac[i, j, k] *
                met[i, j, k + 1, 1, 3] *
                (rhoold[i, j, k + 1] + rhobar[i, j, k + 1]) *
                vu
            ) / (jac[i, j, k] + jac[i, j, k + 1]) -
            fc * (
                jac[i, j, k + 1] *
                met[i, j, k, 2, 3] *
                (rhoold[i, j, k] + rhobar[i, j, k]) *
                uc +
                jac[i, j, k] *
                met[i, j, k + 1, 2, 3] *
                (rhoold[i, j, k + 1] + rhobar[i, j, k + 1]) *
                uu
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        force = -fluxdiff + volforce

        # Interpolate densities.
        rhom_1 =
            (
                jac[i, j, k + 1] * rhoold[i, j, k] +
                jac[i, j, k] * rhoold[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhobaredgeu =
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhom_1 += rhobaredgeu
        rhom += rhobaredgeu

        wm_1 = w[i, j, k]
        momm_1 = rhom_1 * wm_1

        dw[i, j, k] = dt * force + alphark[m] * dw[i, j, k]

        momm = momm_1 + betark[m] * dw[i, j, k]

        # Update wind.
        wast = momm / rhom
        w[i, j, k] = wast
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Explicit,
)
    (; g_ndim) = state.constants
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac) = state.grid
    (; rhobar) = state.atmosphere
    (; rhopold) = state.variables.backups
    (; rho, w, pip) = state.variables.predictands

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 - 1 : k1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        rhoc = rho[i, j, k]
        rhou = rho[i, j, k + 1]
        rhoedgeu =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhoc += rhobar[i, j, k]
        rhou += rhobar[i, j, k + 1]
        rhoedgeu +=
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        gradient = compute_pressure_gradient(state, pip, i, j, k, W())

        force = compute_volume_force(state, i, j, k, W())

        b =
            -g_ndim * (
                jac[i, j, k + 1] * rhopold[i, j, k] / rhoc / jac[i, j, k] +
                jac[i, j, k] * rhopold[i, j, k + 1] / rhou / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jpedgeu = compute_compressible_wind_factor(state, i, j, k, W())

        w[i, j, k] += dt * (b - gradient + force / rhoedgeu) * jpedgeu
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::W,
    side::RHS,
    integration::Implicit,
    rayleigh_factor::AbstractFloat,
)
    (; g_ndim) = state.constants
    (; zz_size, nzz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; jac, met) = state.grid
    (; rhobar, n2) = state.atmosphere
    (; betar) = state.sponge
    (; rho, rhop, u, v, w, pip) = state.variables.predictands

    kmin = ko == 0 ? k0 : k0 - 1
    kmax = ko + nzz == zz_size ? k1 - 1 : k1

    @ivy for k in kmin:kmax, j in j0:j1, i in i0:i1
        rhoc = rho[i, j, k]
        rhou = rho[i, j, k + 1]
        rhoedgeu =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        rhoc += rhobar[i, j, k]
        rhou += rhobar[i, j, k + 1]
        rhoedgeu +=
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        gradient = compute_pressure_gradient(state, pip, i, j, k, W())

        force = compute_volume_force(state, i, j, k, W())

        n2edgeu =
            (jac[i, j, k + 1] * n2[i, j, k] + jac[i, j, k] * n2[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])

        factor = 1.0

        factor +=
            dt * (
                jac[i, j, k + 1] * betar[i, j, k] +
                jac[i, j, k] * betar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1]) * rayleigh_factor

        # Buoyancy is predicted after momentum in implicit steps.
        b =
            -g_ndim * (
                jac[i, j, k + 1] * rhop[i, j, k] / rhoc / jac[i, j, k] +
                jac[i, j, k] * rhop[i, j, k + 1] / rhou / jac[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])

        jpedger = compute_compressible_wind_factor(state, i, j, k, U())
        jpedgel = compute_compressible_wind_factor(state, i - 1, j, k, U())
        jpedgef = compute_compressible_wind_factor(state, i, j, k, V())
        jpedgeb = compute_compressible_wind_factor(state, i, j - 1, k, V())
        jpuedger = compute_compressible_wind_factor(state, i, j, k + 1, U())
        jpuedgel = compute_compressible_wind_factor(state, i - 1, j, k + 1, U())
        jpuedgef = compute_compressible_wind_factor(state, i, j, k + 1, V())
        jpuedgeb = compute_compressible_wind_factor(state, i, j - 1, k + 1, V())

        uc = 0.5 * (u[i, j, k] / jpedger + u[i - 1, j, k] / jpedgel)
        uu = 0.5 * (u[i, j, k + 1] / jpuedger + u[i - 1, j, k + 1] / jpuedgel)
        vc = 0.5 * (v[i, j, k] / jpedgef + v[i, j - 1, k] / jpedgeb)
        vu = 0.5 * (v[i, j, k + 1] / jpuedgef + v[i, j - 1, k + 1] / jpuedgeb)

        jpedgeu = compute_compressible_wind_factor(state, i, j, k, W())
        fw = compute_buoyancy_factor(state, i, j, k, W())

        w[i, j, k] =
            1.0 / (factor + fw * n2edgeu * dt^2.0) * (
                w[i, j, k] - dt * gradient * jpedgeu +
                dt * b * jpedgeu +
                dt * force / rhoedgeu * jpedgeu +
                jpedgeu *
                fw *
                n2edgeu *
                dt^2.0 *
                (
                    jac[i, j, k + 1] *
                    (met[i, j, k, 1, 3] * uc + met[i, j, k, 2, 3] * vc) +
                    jac[i, j, k] *
                    (met[i, j, k + 1, 1, 3] * uu + met[i, j, k + 1, 2, 3] * vu)
                ) / (jac[i, j, k] + jac[i, j, k + 1])
            )
    end

    return
end

function update!(state::State, dt::AbstractFloat, variable::PiP)
    (; model) = state.namelists.setting
    update!(state, dt, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::PiP,
    model::Union{Boussinesq, PseudoIncompressible},
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    variable::PiP,
    model::Compressible,
)
    (; gamma, rsp, pref) = state.constants
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; uold, vold, wold) = state.variables.backups
    (; pip, p) = state.variables.predictands

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        fl = uold[i - 1, j, k]
        fr = uold[i, j, k]
        gb = vold[i, j - 1, k]
        gf = vold[i, j, k]
        hd = wold[i, j, k - 1]
        hu = wold[i, j, k]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, i, j, k, P())

        dpdpi =
            1 / (gamma - 1) * (rsp / pref)^(1 - gamma) * p[i, j, k]^(2 - gamma)

        pip[i, j, k] -= dt * (fluxdiff - heating) / dpdpi
    end

    return
end

function update!(state::State, dt::AbstractFloat, m::Integer, variable::P)
    (; model) = state.namelists.setting
    update!(state, dt, m, variable, model)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Union{Boussinesq, PseudoIncompressible},
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    variable::P,
    model::Compressible,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; dp) = state.variables.increments
    (; phip) = state.variables.fluxes
    (; p) = state.variables.predictands

    if m == 1
        dp .= 0.0
    end

    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        fl = phip[i - 1, j, k, 1]
        fr = phip[i, j, k, 1]
        gb = phip[i, j - 1, k, 2]
        gf = phip[i, j, k, 2]
        hd = phip[i, j, k - 1, 3]
        hu = phip[i, j, k, 3]

        fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
        fluxdiff /= jac[i, j, k]

        heating = compute_volume_force(state, i, j, k, P())

        f = -fluxdiff + heating

        dp[i, j, k] = dt * f + alphark[m] * dp[i, j, k]
        p[i, j, k] += betark[m] * dp[i, j, k]
    end

    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    tracer_setup::NoTracer,
)
    return
end

function update!(
    state::State,
    dt::AbstractFloat,
    m::Integer,
    tracer_setup::TracerOn,
)
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac) = state.grid
    (; alphark, betark) = state.time
    (; tracerincrements, tracerpredictands, tracerfluxes) = state.tracer

    @ivy for field in 1:fieldcount(TracerPredictands)
        if m == 1
            getfield(tracerincrements, field) .= 0.0
        end

        flr = getfield(tracerfluxes, field)[:, :, :, 1]
        gbf = getfield(tracerfluxes, field)[:, :, :, 2]
        hdu = getfield(tracerfluxes, field)[:, :, :, 3]
        chi = getfield(tracerpredictands, field)[:, :, :]
        dchi = getfield(tracerincrements, field)[:, :, :]
        for k in k0:k1, j in j0:j1, i in i0:i1
            fl = flr[i - 1, j, k]
            fr = flr[i, j, k]
            gb = gbf[i, j - 1, k]
            gf = gbf[i, j, k]
            hd = hdu[i, j, k - 1]
            hu = hdu[i, j, k]

            fluxdiff = (fr - fl) / dx + (gf - gb) / dy + (hu - hd) / dz
            fluxdiff /= jac[i, j, k]

            force = compute_volume_force(state, i, j, k, Chi())
            f = -fluxdiff + force

            dchi[i, j, k] = dt * f + alphark[m] * dchi[i, j, k]
            chi[i, j, k] += betark[m] * dchi[i, j, k]
        end
    end

    return
end
