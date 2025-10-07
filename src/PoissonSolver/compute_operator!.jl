"""
```julia
compute_operator!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)
```

Compute the tensor elements of the linear operator on the right-hand side of the Poisson equation.

The operator is obtained by rewriting the scaled Poisson equation

```math
\\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{LHS} = \\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{RHS} \\left(\\frac{\\sqrt{\\overline{\\rho}}}{P} s\\right)
```

as

```math
\\frac{\\sqrt{\\overline{\\rho}}}{P} \\mathrm{LHS} = \\sum_{\\lambda, \\mu, \\nu} A_{i + \\lambda, j + \\mu, k + \\nu} s_{i + \\lambda, j + \\mu, k + \\nu},
```

where the Exner-pressure differences are given by ``\\Delta \\pi = \\left(\\sqrt{\\overline{\\rho}} / P\\right) \\left(s / \\Delta t\\right)``.

# Arguments

  - `state`: Model state.

  - `dt`: Time step.

  - `rayleigh_factor`: Factor by which the Rayleigh-damping coefficient is multiplied.

# See also

  - [`PinCFlow.Update.compute_buoyancy_factor`](@ref)
"""
function compute_operator! end

function compute_operator!(
    state::State,
    dt::AbstractFloat,
    rayleigh_factor::AbstractFloat,
)
    (; nbz) = state.namelists.domain
    (; spongelayer, sponge_uv) = state.namelists.sponge
    (; model) = state.namelists.setting
    (; gamma, rsp, pref, kappa) = state.constants
    (; sizezz, ko, i0, i1, j0, j1, k0, k1) = state.domain
    (; dx, dy, dz, jac, met) = state.grid
    (; pbar, rhobar, n2) = state.atmosphere
    (; betar) = state.sponge
    (;
        ac_b,
        al_b,
        ar_b,
        ab_b,
        af_b,
        ad_b,
        au_b,
        aru_b,
        ard_b,
        alu_b,
        ald_b,
        afu_b,
        afd_b,
        abu_b,
        abd_b,
        auu_b,
        add_b,
        aruu_b,
        ardd_b,
        aluu_b,
        aldd_b,
        afuu_b,
        afdd_b,
        abuu_b,
        abdd_b,
    ) = state.poisson.tensor
    (; rho, p) = state.variables.predictands

    # Compute tensor elements for TFC.
    @ivy for k in k0:k1, j in j0:j1, i in i0:i1
        # Compute scaling factors.
        fcscal = sqrt(pbar[i, j, k]^2.0 / rhobar[i, j, k])
        fcscal_r = sqrt(pbar[i + 1, j, k]^2.0 / rhobar[i + 1, j, k])
        fcscal_l = sqrt(pbar[i - 1, j, k]^2.0 / rhobar[i - 1, j, k])
        fcscal_f = sqrt(pbar[i, j + 1, k]^2.0 / rhobar[i, j + 1, k])
        fcscal_b = sqrt(pbar[i, j - 1, k]^2.0 / rhobar[i, j - 1, k])
        fcscal_u = sqrt(pbar[i, j, k + 1]^2.0 / rhobar[i, j, k + 1])
        fcscal_d = sqrt(pbar[i, j, k - 1]^2.0 / rhobar[i, j, k - 1])
        fcscal_ru = sqrt(pbar[i + 1, j, k + 1]^2.0 / rhobar[i + 1, j, k + 1])
        fcscal_rd = sqrt(pbar[i + 1, j, k - 1]^2.0 / rhobar[i + 1, j, k - 1])
        fcscal_lu = sqrt(pbar[i - 1, j, k + 1]^2.0 / rhobar[i - 1, j, k + 1])
        fcscal_ld = sqrt(pbar[i - 1, j, k - 1]^2.0 / rhobar[i - 1, j, k - 1])
        fcscal_fu = sqrt(pbar[i, j + 1, k + 1]^2.0 / rhobar[i, j + 1, k + 1])
        fcscal_fd = sqrt(pbar[i, j + 1, k - 1]^2.0 / rhobar[i, j + 1, k - 1])
        fcscal_bu = sqrt(pbar[i, j - 1, k + 1]^2.0 / rhobar[i, j - 1, k + 1])
        fcscal_bd = sqrt(pbar[i, j - 1, k - 1]^2.0 / rhobar[i, j - 1, k - 1])
        fcscal_uu = sqrt(pbar[i, j, k + 2]^2.0 / rhobar[i, j, k + 2])
        fcscal_dd = sqrt(pbar[i, j, k - 2]^2.0 / rhobar[i, j, k - 2])
        fcscal_ruu = sqrt(pbar[i + 1, j, k + 2]^2.0 / rhobar[i + 1, j, k + 2])
        fcscal_rdd = sqrt(pbar[i + 1, j, k - 2]^2.0 / rhobar[i + 1, j, k - 2])
        fcscal_luu = sqrt(pbar[i - 1, j, k + 2]^2.0 / rhobar[i - 1, j, k + 2])
        fcscal_ldd = sqrt(pbar[i - 1, j, k - 2]^2.0 / rhobar[i - 1, j, k - 2])
        fcscal_fuu = sqrt(pbar[i, j + 1, k + 2]^2.0 / rhobar[i, j + 1, k + 2])
        fcscal_fdd = sqrt(pbar[i, j + 1, k - 2]^2.0 / rhobar[i, j + 1, k - 2])
        fcscal_buu = sqrt(pbar[i, j - 1, k + 2]^2.0 / rhobar[i, j - 1, k + 2])
        fcscal_bdd = sqrt(pbar[i, j - 1, k - 2]^2.0 / rhobar[i, j - 1, k - 2])

        # Compute inverse Jacobian.
        jacinv = 1.0 / jac[i, j, k]

        # Compute P coefficients (divergence).
        pedgerdiv =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i + 1, j, k] * pbar[i + 1, j, k]
            )
        pedgeldiv =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i - 1, j, k] * pbar[i - 1, j, k]
            )
        pedgefdiv =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j + 1, k] * pbar[i, j + 1, k]
            )
        pedgebdiv =
            0.5 * (
                jac[i, j, k] * pbar[i, j, k] +
                jac[i, j - 1, k] * pbar[i, j - 1, k]
            )
        pedgeudiv =
            jac[i, j, k] *
            jac[i, j, k + 1] *
            (pbar[i, j, k] + pbar[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        pedgeddiv =
            jac[i, j, k] *
            jac[i, j, k - 1] *
            (pbar[i, j, k] + pbar[i, j, k - 1]) /
            (jac[i, j, k] + jac[i, j, k - 1])

        # Compute P coefficients (pressure gradient).
        pedgergra = 0.5 * (pbar[i, j, k] + pbar[i + 1, j, k])
        pedgelgra = 0.5 * (pbar[i, j, k] + pbar[i - 1, j, k])
        pedgefgra = 0.5 * (pbar[i, j, k] + pbar[i, j + 1, k])
        pedgebgra = 0.5 * (pbar[i, j, k] + pbar[i, j - 1, k])
        pedgeugra =
            (
                jac[i, j, k + 1] * pbar[i, j, k] +
                jac[i, j, k] * pbar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        pedgedgra =
            (
                jac[i, j, k - 1] * pbar[i, j, k] +
                jac[i, j, k] * pbar[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        puedgergra = 0.5 * (pbar[i, j, k + 1] + pbar[i + 1, j, k + 1])
        puedgelgra = 0.5 * (pbar[i, j, k + 1] + pbar[i - 1, j, k + 1])
        puedgefgra = 0.5 * (pbar[i, j, k + 1] + pbar[i, j + 1, k + 1])
        puedgebgra = 0.5 * (pbar[i, j, k + 1] + pbar[i, j - 1, k + 1])
        pdedgergra = 0.5 * (pbar[i, j, k - 1] + pbar[i + 1, j, k - 1])
        pdedgelgra = 0.5 * (pbar[i, j, k - 1] + pbar[i - 1, j, k - 1])
        pdedgefgra = 0.5 * (pbar[i, j, k - 1] + pbar[i, j + 1, k - 1])
        pdedgebgra = 0.5 * (pbar[i, j, k - 1] + pbar[i, j - 1, k - 1])

        # Compute density coefficients.
        rhobaredger = 0.5 * (rhobar[i, j, k] + rhobar[i + 1, j, k])
        rhobaredgel = 0.5 * (rhobar[i, j, k] + rhobar[i - 1, j, k])
        rhobaredgef = 0.5 * (rhobar[i, j, k] + rhobar[i, j + 1, k])
        rhobaredgeb = 0.5 * (rhobar[i, j, k] + rhobar[i, j - 1, k])
        rhobaredgeu =
            (
                jac[i, j, k + 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        rhobaredged =
            (
                jac[i, j, k - 1] * rhobar[i, j, k] +
                jac[i, j, k] * rhobar[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        rhoedger = 0.5 * (rho[i, j, k] + rho[i + 1, j, k]) + rhobaredger
        rhoedgel = 0.5 * (rho[i, j, k] + rho[i - 1, j, k]) + rhobaredgel
        rhoedgef = 0.5 * (rho[i, j, k] + rho[i, j + 1, k]) + rhobaredgef
        rhoedgeb = 0.5 * (rho[i, j, k] + rho[i, j - 1, k]) + rhobaredgeb
        rhoedgeu =
            (
                jac[i, j, k + 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k + 1]
            ) / (jac[i, j, k] + jac[i, j, k + 1]) + rhobaredgeu
        rhoedged =
            (
                jac[i, j, k - 1] * rho[i, j, k] +
                jac[i, j, k] * rho[i, j, k - 1]
            ) / (jac[i, j, k] + jac[i, j, k - 1]) + rhobaredged

        fwu = compute_buoyancy_factor(state, i, j, k, W())
        fwd = compute_buoyancy_factor(state, i, j, k - 1, W())

        rhouedger =
            0.5 * (
                rho[i, j, k + 1] +
                rho[i + 1, j, k + 1] +
                rhobar[i, j, k + 1] +
                rhobar[i + 1, j, k + 1]
            )
        rhouedgel =
            0.5 * (
                rho[i, j, k + 1] +
                rho[i - 1, j, k + 1] +
                rhobar[i, j, k + 1] +
                rhobar[i - 1, j, k + 1]
            )
        rhouedgef =
            0.5 * (
                rho[i, j, k + 1] +
                rho[i, j + 1, k + 1] +
                rhobar[i, j, k + 1] +
                rhobar[i, j + 1, k + 1]
            )
        rhouedgeb =
            0.5 * (
                rho[i, j, k + 1] +
                rho[i, j - 1, k + 1] +
                rhobar[i, j, k + 1] +
                rhobar[i, j - 1, k + 1]
            )
        rhodedger =
            0.5 * (
                rho[i, j, k - 1] +
                rho[i + 1, j, k - 1] +
                rhobar[i, j, k - 1] +
                rhobar[i + 1, j, k - 1]
            )
        rhodedgel =
            0.5 * (
                rho[i, j, k - 1] +
                rho[i - 1, j, k - 1] +
                rhobar[i, j, k - 1] +
                rhobar[i - 1, j, k - 1]
            )
        rhodedgef =
            0.5 * (
                rho[i, j, k - 1] +
                rho[i, j + 1, k - 1] +
                rhobar[i, j, k - 1] +
                rhobar[i, j + 1, k - 1]
            )
        rhodedgeb =
            0.5 * (
                rho[i, j, k - 1] +
                rho[i, j - 1, k - 1] +
                rhobar[i, j, k - 1] +
                rhobar[i, j - 1, k - 1]
            )

        # Compute squared buoyancy frequency at edges.
        n2edgeu =
            (jac[i, j, k + 1] * n2[i, j, k] + jac[i, j, k] * n2[i, j, k + 1]) /
            (jac[i, j, k] + jac[i, j, k + 1])
        n2edged =
            (jac[i, j, k - 1] * n2[i, j, k] + jac[i, j, k] * n2[i, j, k - 1]) /
            (jac[i, j, k] + jac[i, j, k - 1])

        # Interpolate metric-tensor elements.
        met13edger = 0.5 * (met[i, j, k, 1, 3] + met[i + 1, j, k, 1, 3])
        met13edgel = 0.5 * (met[i, j, k, 1, 3] + met[i - 1, j, k, 1, 3])
        met23edgef = 0.5 * (met[i, j, k, 2, 3] + met[i, j + 1, k, 2, 3])
        met23edgeb = 0.5 * (met[i, j, k, 2, 3] + met[i, j - 1, k, 2, 3])
        met13edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k + 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met23edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k + 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met33edgeu =
            (
                jac[i, j, k + 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k + 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k + 1])
        met13edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 1, 3] +
                jac[i, j, k] * met[i, j, k - 1, 1, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met23edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 2, 3] +
                jac[i, j, k] * met[i, j, k - 1, 2, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met33edged =
            (
                jac[i, j, k - 1] * met[i, j, k, 3, 3] +
                jac[i, j, k] * met[i, j, k - 1, 3, 3]
            ) / (jac[i, j, k] + jac[i, j, k - 1])
        met13uedger =
            0.5 * (met[i, j, k + 1, 1, 3] + met[i + 1, j, k + 1, 1, 3])
        met13uedgel =
            0.5 * (met[i, j, k + 1, 1, 3] + met[i - 1, j, k + 1, 1, 3])
        met23uedgef =
            0.5 * (met[i, j, k + 1, 2, 3] + met[i, j + 1, k + 1, 2, 3])
        met23uedgeb =
            0.5 * (met[i, j, k + 1, 2, 3] + met[i, j - 1, k + 1, 2, 3])
        met13dedger =
            0.5 * (met[i, j, k - 1, 1, 3] + met[i + 1, j, k - 1, 1, 3])
        met13dedgel =
            0.5 * (met[i, j, k - 1, 1, 3] + met[i - 1, j, k - 1, 1, 3])
        met23dedgef =
            0.5 * (met[i, j, k - 1, 2, 3] + met[i, j + 1, k - 1, 2, 3])
        met23dedgeb =
            0.5 * (met[i, j, k - 1, 2, 3] + met[i, j - 1, k - 1, 2, 3])

        # Compute Rayleigh damping terms.
        facedger = 1.0
        facedgel = 1.0
        facedgef = 1.0
        facedgeb = 1.0
        facuedger = 1.0
        facuedgel = 1.0
        facuedgef = 1.0
        facuedgeb = 1.0
        facdedger = 1.0
        facdedgel = 1.0
        facdedgef = 1.0
        facdedgeb = 1.0
        facedgeu = 1.0
        facedged = 1.0
        if spongelayer
            if sponge_uv
                facedger =
                    facedger +
                    dt *
                    0.5 *
                    (betar[i, j, k] + betar[i + 1, j, k]) *
                    rayleigh_factor
                facedgel =
                    facedgel +
                    dt *
                    0.5 *
                    (betar[i, j, k] + betar[i - 1, j, k]) *
                    rayleigh_factor
                facedgef =
                    facedgef +
                    dt *
                    0.5 *
                    (betar[i, j, k] + betar[i, j + 1, k]) *
                    rayleigh_factor
                facedgeb =
                    facedgeb +
                    dt *
                    0.5 *
                    (betar[i, j, k] + betar[i, j - 1, k]) *
                    rayleigh_factor
                facuedger =
                    facuedger +
                    dt *
                    0.5 *
                    (betar[i, j, k + 1] + betar[i + 1, j, k + 1]) *
                    rayleigh_factor
                facuedgel =
                    facuedgel +
                    dt *
                    0.5 *
                    (betar[i, j, k + 1] + betar[i - 1, j, k + 1]) *
                    rayleigh_factor
                facuedgef =
                    facuedgef +
                    dt *
                    0.5 *
                    (betar[i, j, k + 1] + betar[i, j + 1, k + 1]) *
                    rayleigh_factor
                facuedgeb =
                    facuedgeb +
                    dt *
                    0.5 *
                    (betar[i, j, k + 1] + betar[i, j - 1, k + 1]) *
                    rayleigh_factor
                facdedger =
                    facdedger +
                    dt *
                    0.5 *
                    (betar[i, j, k - 1] + betar[i + 1, j, k - 1]) *
                    rayleigh_factor
                facdedgel =
                    facdedgel +
                    dt *
                    0.5 *
                    (betar[i, j, k - 1] + betar[i - 1, j, k - 1]) *
                    rayleigh_factor
                facdedgef =
                    facdedgef +
                    dt *
                    0.5 *
                    (betar[i, j, k - 1] + betar[i, j + 1, k - 1]) *
                    rayleigh_factor
                facdedgeb =
                    facdedgeb +
                    dt *
                    0.5 *
                    (betar[i, j, k - 1] + betar[i, j - 1, k - 1]) *
                    rayleigh_factor
            end
            facedgeu =
                facedgeu +
                dt * (
                    jac[i, j, k + 1] * betar[i, j, k] +
                    jac[i, j, k] * betar[i, j, k + 1]
                ) / (jac[i, j, k] + jac[i, j, k + 1]) * rayleigh_factor
            facedged =
                facedged +
                dt * (
                    jac[i, j, k - 1] * betar[i, j, k] +
                    jac[i, j, k] * betar[i, j, k - 1]
                ) / (jac[i, j, k] + jac[i, j, k - 1]) * rayleigh_factor
        end

        # Compute implicit coefficients.
        imphoredger = 1.0 / facedger
        imphoredgel = 1.0 / facedgel
        imphoredgef = 1.0 / facedgef
        imphoredgeb = 1.0 / facedgeb
        imphoruedger = 1.0 / facuedger
        imphoruedgel = 1.0 / facuedgel
        imphoruedgef = 1.0 / facuedgef
        imphoruedgeb = 1.0 / facuedgeb
        imphordedger = 1.0 / facdedger
        imphordedgel = 1.0 / facdedgel
        imphordedgef = 1.0 / facdedgef
        imphordedgeb = 1.0 / facdedgeb
        impveredgeu = 1.0 / (facedgeu + fwu * n2edgeu * dt^2.0)
        impveredged = 1.0 / (facedged + fwd * n2edged * dt^2.0)

        # Compute gradient coefficients

        # G(i + 1 / 2)
        if ko + k == k0
            gedger =
                jacinv / dx * pedgerdiv * imphoredger / rhoedger +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredger / rhoedger
        elseif ko + k == sizezz - nbz
            gedger =
                jacinv / dx * pedgerdiv * imphoredger / rhoedger -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredger / rhoedger
        else
            gedger =
                jacinv / dx * pedgerdiv * imphoredger / rhoedger +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredger / rhoedger -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredger / rhoedger
        end

        # G(i - 1 / 2)
        if ko + k == k0
            gedgel =
                -jacinv / dx * pedgeldiv * imphoredgel / rhoedgel +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgel / rhoedgel
        elseif ko + k == sizezz - nbz
            gedgel =
                -jacinv / dx * pedgeldiv * imphoredgel / rhoedgel -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgel / rhoedgel
        else
            gedgel =
                -jacinv / dx * pedgeldiv * imphoredgel / rhoedgel +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgel / rhoedgel -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 1, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgel / rhoedgel
        end

        # G(j + 1 / 2)
        if ko + k == k0
            gedgef =
                jacinv / dy * pedgefdiv * imphoredgef / rhoedgef +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgef / rhoedgef
        elseif ko + k == sizezz - nbz
            gedgef =
                jacinv / dy * pedgefdiv * imphoredgef / rhoedgef -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgef / rhoedgef
        else
            gedgef =
                jacinv / dy * pedgefdiv * imphoredgef / rhoedgef +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgef / rhoedgef -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgef / rhoedgef
        end

        # G(j - 1 / 2)
        if ko + k == k0
            gedgeb =
                -jacinv / dy * pedgebdiv * imphoredgeb / rhoedgeb +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgeb / rhoedgeb
        elseif ko + k == sizezz - nbz
            gedgeb =
                -jacinv / dy * pedgebdiv * imphoredgeb / rhoedgeb -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgeb / rhoedgeb
        else
            gedgeb =
                -jacinv / dy * pedgebdiv * imphoredgeb / rhoedgeb +
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k + 1] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoredgeb / rhoedgeb -
                jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k, 2, 3] *
                jac[i, j, k - 1] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphoredgeb / rhoedgeb
        end

        # G(k + 1 / 2)
        if ko + k == sizezz - nbz
            gedgeu = 0.0
        else
            gedgeu = jacinv / dz * pedgeudiv * impveredgeu / rhoedgeu
        end

        # G(k - 1 / 2)
        if ko + k == k0
            gedged = 0.0
        else
            gedged = -jacinv / dz * pedgeddiv * impveredged / rhoedged
        end

        # G(i + 1 / 2, k + 1)
        if ko + k == sizezz - nbz
            guedger = 0.0
        else
            guedger =
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k + 1, 1, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoruedger / rhouedger
        end

        # G(i - 1 / 2, k + 1)
        if ko + k == sizezz - nbz
            guedgel = 0.0
        else
            guedgel =
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k + 1, 1, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoruedgel / rhouedgel
        end

        # G(j + 1 / 2, k + 1)
        if ko + k == sizezz - nbz
            guedgef = 0.0
        else
            guedgef =
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k + 1, 2, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoruedgef / rhouedgef
        end

        # G(j - 1 / 2, k + 1)
        if ko + k == sizezz - nbz
            guedgeb = 0.0
        else
            guedgeb =
                jacinv / dz *
                pedgeudiv *
                impveredgeu *
                fwu *
                n2edgeu *
                dt^2.0 *
                0.5 *
                met[i, j, k + 1, 2, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k + 1]) *
                imphoruedgeb / rhouedgeb
        end

        # G(i + 1 / 2, k - 1)
        if ko + k == k0
            gdedger = 0.0
        else
            gdedger =
                -jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k - 1, 1, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphordedger / rhodedger
        end

        # G(i - 1 / 2, k - 1)
        if ko + k == k0
            gdedgel = 0.0
        else
            gdedgel =
                -jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k - 1, 1, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphordedgel / rhodedgel
        end

        # G(j + 1 / 2, k - 1)
        if ko + k == k0
            gdedgef = 0.0
        else
            gdedgef =
                -jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k - 1, 2, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphordedgef / rhodedgef
        end

        # G(j - 1 / 2, k - 1)
        if ko + k == k0
            gdedgeb = 0.0
        else
            gdedgeb =
                -jacinv / dz *
                pedgeddiv *
                impveredged *
                fwd *
                n2edged *
                dt^2.0 *
                0.5 *
                met[i, j, k - 1, 2, 3] *
                jac[i, j, k] / (jac[i, j, k] + jac[i, j, k - 1]) *
                imphordedgeb / rhodedgeb
        end

        # Compute tensor elements

        # ------------------- A(i,j,k) --------------------#

        if ko + k == k0
            ac =
                -gedger * pedgergra * (1.0 / dx + 0.75 * met13edger / dz) +
                gedgel * pedgelgra * (1.0 / dx - 0.75 * met13edgel / dz) -
                gedgef * pedgefgra * (1.0 / dy + 0.75 * met23edgef / dz) +
                gedgeb * pedgebgra * (1.0 / dy - 0.75 * met23edgeb / dz) -
                gedgeu * pedgeugra * met33edgeu / dz -
                guedger * puedgergra * 0.25 * met13uedger / dz -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k == k0 + 1
            ac =
                -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
                gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * met33edgeu / dz +
                gedged * pedgedgra * met33edged / dz -
                guedger * puedgergra * 0.25 * met13uedger / dz -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
                gdedger * pdedgergra * met13dedger / dz +
                gdedgel * pdedgelgra * met13dedgel / dz +
                gdedgef * pdedgefgra * met23dedgef / dz +
                gdedgeb * pdedgebgra * met23dedgeb / dz
        elseif ko + k == sizezz - nbz - 1
            ac =
                -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
                gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * met33edgeu / dz +
                gedged * pedgedgra * met33edged / dz -
                guedger * puedgergra * met13uedger / dz -
                guedgel * puedgelgra * met13uedgel / dz -
                guedgef * puedgefgra * met23uedgef / dz -
                guedgeb * puedgebgra * met23uedgeb / dz +
                gdedger * pdedgergra * 0.25 * met13dedger / dz +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        elseif ko + k == sizezz - nbz
            ac =
                -gedger * pedgergra * (1.0 / dx - 0.75 * met13edger / dz) +
                gedgel * pedgelgra * (1.0 / dx + 0.75 * met13edgel / dz) -
                gedgef * pedgefgra * (1.0 / dy - 0.75 * met23edgef / dz) +
                gedgeb * pedgebgra * (1.0 / dy + 0.75 * met23edgeb / dz) +
                gedged * pedgedgra * met33edged / dz +
                gdedger * pdedgergra * 0.25 * met13dedger / dz +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        else
            ac =
                -gedger * pedgergra / dx + gedgel * pedgelgra / dx -
                gedgef * pedgefgra / dy + gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * met33edgeu / dz +
                gedged * pedgedgra * met33edged / dz -
                guedger * puedgergra * 0.25 * met13uedger / dz -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
                gdedger * pdedgergra * 0.25 * met13dedger / dz +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        end

        if model == Compressible()
            dpdpi =
                1 / (gamma - 1) *
                (rsp / pref)^(1 - gamma) *
                p[i, j, k]^(2 - gamma)
            ac -= (dpdpi / dt) / (dt * rsp / kappa)
        end

        # ------------------ A(i+1,j,k) -------------------#

        if ko + k == k0
            ar =
                gedger * pedgergra * (1.0 / dx - 0.75 * met13edger / dz) +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i + 1, j, k + 1] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) -
                guedger * puedgergra * 0.25 * met13uedger / dz
        elseif ko + k == k0 + 1
            ar =
                gedger * pedgergra / dx +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i + 1, j, k + 1] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i + 1, j, k - 1] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
                guedger * puedgergra * 0.25 * met13uedger / dz +
                gdedger * pdedgergra * met13dedger / dz
        elseif ko + k == sizezz - nbz - 1
            ar =
                gedger * pedgergra / dx +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i + 1, j, k + 1] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i + 1, j, k - 1] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
                guedger * puedgergra * met13uedger / dz +
                gdedger * pdedgergra * 0.25 * met13dedger / dz
        elseif ko + k == sizezz - nbz
            ar =
                gedger * pedgergra * (1.0 / dx + 0.75 * met13edger / dz) +
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i + 1, j, k - 1] /
                (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
                gdedger * pdedgergra * 0.25 * met13dedger / dz
        else
            ar =
                gedger * pedgergra / dx +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i + 1, j, k + 1] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i + 1, j, k - 1] / (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) -
                guedger * puedgergra * 0.25 * met13uedger / dz +
                gdedger * pdedgergra * 0.25 * met13dedger / dz
        end

        # ------------------ A(i-1,j,k) -------------------#

        if ko + k == k0
            al =
                -gedgel * pedgelgra * (1.0 / dx + 0.75 * met13edgel / dz) -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i - 1, j, k + 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz
        elseif ko + k == k0 + 1
            al =
                -gedgel * pedgelgra / dx -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i - 1, j, k + 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i - 1, j, k - 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz +
                gdedgel * pdedgelgra * met13dedgel / dz
        elseif ko + k == sizezz - nbz - 1
            al =
                -gedgel * pedgelgra / dx -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i - 1, j, k + 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i - 1, j, k - 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                guedgel * puedgelgra * met13uedgel / dz +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        elseif ko + k == sizezz - nbz
            al =
                -gedgel * pedgelgra * (1.0 / dx - 0.75 * met13edgel / dz) -
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i - 1, j, k - 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        else
            al =
                -gedgel * pedgelgra / dx -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx *
                jac[i - 1, j, k + 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                gedged * pedgedgra * 0.5 * met13edged / dx *
                jac[i - 1, j, k - 1] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                guedgel * puedgelgra * 0.25 * met13uedgel / dz +
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        end

        # ------------------ A(i,j+1,k) -------------------#

        if ko + k == k0
            af =
                gedgef * pedgefgra * (1.0 / dy - 0.75 * met23edgef / dz) +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j + 1, k + 1] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz
        elseif ko + k == k0 + 1
            af =
                gedgef * pedgefgra / dy +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j + 1, k + 1] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j + 1, k - 1] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz +
                gdedgef * pdedgefgra * met23dedgef / dz
        elseif ko + k == sizezz - nbz - 1
            af =
                gedgef * pedgefgra / dy +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j + 1, k + 1] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j + 1, k - 1] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
                guedgef * puedgefgra * met23uedgef / dz +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        elseif ko + k == sizezz - nbz
            af =
                gedgef * pedgefgra * (1.0 / dy + 0.75 * met23edgef / dz) +
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j + 1, k - 1] /
                (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        else
            af =
                gedgef * pedgefgra / dy +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j + 1, k + 1] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j + 1, k - 1] / (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) -
                guedgef * puedgefgra * 0.25 * met23uedgef / dz +
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        end

        # ------------------ A(i,j-1,k) -------------------#

        if ko + k == k0
            ab =
                -gedgeb * pedgebgra * (1.0 / dy + 0.75 * met23edgeb / dz) -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j - 1, k + 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k == k0 + 1
            ab =
                -gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j - 1, k + 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j - 1, k - 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
                gdedgeb * pdedgebgra * met23dedgeb / dz
        elseif ko + k == sizezz - nbz - 1
            ab =
                -gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j - 1, k + 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j - 1, k - 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                guedgeb * puedgebgra * met23uedgeb / dz +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        elseif ko + k == sizezz - nbz
            ab =
                -gedgeb * pedgebgra * (1.0 / dy - 0.75 * met23edgeb / dz) -
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j - 1, k - 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        else
            ab =
                -gedgeb * pedgebgra / dy -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy *
                jac[i, j - 1, k + 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                gedged * pedgedgra * 0.5 * met23edged / dy *
                jac[i, j - 1, k - 1] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz +
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        end

        # ------------------ A(i,j,k+1) -------------------#

        if ko + k == k0
            au =
                gedger * pedgergra * met13edger / dz +
                gedgel * pedgelgra * met13edgel / dz +
                gedgef * pedgefgra * met23edgef / dz +
                gedgeb * pedgebgra * met23edgeb / dz +
                gedgeu * pedgeugra * met33edgeu / dz -
                guedger * puedgergra / dx + guedgel * puedgelgra / dx -
                guedgef * puedgefgra / dy + guedgeb * puedgebgra / dy
        elseif ko + k == k0 + 1
            au =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgel * pedgelgra * 0.25 * met13edgel / dz +
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
                gedgeu * pedgeugra * met33edgeu / dz -
                guedger * puedgergra / dx + guedgel * puedgelgra / dx -
                guedgef * puedgefgra / dy + guedgeb * puedgebgra / dy -
                gdedger * pdedgergra * 0.25 * met13dedger / dz -
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        elseif ko + k == sizezz - nbz - 1
            au =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgel * pedgelgra * 0.25 * met13edgel / dz +
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
                gedgeu * pedgeugra * met33edgeu / dz -
                guedger * puedgergra * (1.0 / dx - 0.75 * met13uedger / dz) +
                guedgel * puedgelgra * (1.0 / dx + 0.75 * met13uedgel / dz) -
                guedgef * puedgefgra * (1.0 / dy - 0.75 * met23uedgef / dz) +
                guedgeb * puedgebgra * (1.0 / dy + 0.75 * met23uedgeb / dz)
        elseif ko + k == sizezz - nbz
            au = 0.0
        else
            au =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgel * pedgelgra * 0.25 * met13edgel / dz +
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
                gedgeu * pedgeugra * met33edgeu / dz -
                guedger * puedgergra / dx + guedgel * puedgelgra / dx -
                guedgef * puedgefgra / dy + guedgeb * puedgebgra / dy
        end

        # ------------------ A(i,j,k-1) -------------------#

        if ko + k == k0
            ad = 0.0
        elseif ko + k == k0 + 1
            ad =
                -gedger * pedgergra * 0.25 * met13edger / dz -
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgef * pedgefgra * 0.25 * met23edgef / dz -
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * met33edged / dz -
                gdedger * pdedgergra * (1.0 / dx + 0.75 * met13dedger / dz) +
                gdedgel * pdedgelgra * (1.0 / dx - 0.75 * met13dedgel / dz) -
                gdedgef * pdedgefgra * (1.0 / dy + 0.75 * met23dedgef / dz) +
                gdedgeb * pdedgebgra * (1.0 / dy - 0.75 * met23dedgeb / dz)
        elseif ko + k == sizezz - nbz - 1
            ad =
                -gedger * pedgergra * 0.25 * met13edger / dz -
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgef * pedgefgra * 0.25 * met23edgef / dz -
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * met33edged / dz -
                gdedger * pdedgergra / dx + gdedgel * pdedgelgra / dx -
                gdedgef * pdedgefgra / dy +
                gdedgeb * pdedgebgra / dy +
                guedger * puedgergra * 0.25 * met13uedger / dz +
                guedgel * puedgelgra * 0.25 * met13uedgel / dz +
                guedgef * puedgefgra * 0.25 * met23uedgef / dz +
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k == sizezz - nbz
            ad =
                -gedger * pedgergra * met13edger / dz -
                gedgel * pedgelgra * met13edgel / dz -
                gedgef * pedgefgra * met23edgef / dz -
                gedgeb * pedgebgra * met23edgeb / dz -
                gedged * pedgedgra * met33edged / dz -
                gdedger * pdedgergra / dx + gdedgel * pdedgelgra / dx -
                gdedgef * pdedgefgra / dy + gdedgeb * pdedgebgra / dy
        else
            ad =
                -gedger * pedgergra * 0.25 * met13edger / dz -
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgef * pedgefgra * 0.25 * met23edgef / dz -
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * met33edged / dz -
                gdedger * pdedgergra / dx + gdedgel * pdedgelgra / dx -
                gdedgef * pdedgefgra / dy + gdedgeb * pdedgebgra / dy
        end

        # ----------------- A(i+1,j,k+1) ------------------#

        if ko + k == k0
            aru =
                gedger * pedgergra * met13edger / dz +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                guedger * puedgergra / dx
        elseif ko + k == k0 + 1
            aru =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                guedger * puedgergra / dx -
                gdedger * pdedgergra * 0.25 * met13dedger / dz
        elseif ko + k == sizezz - nbz - 1
            aru =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                guedger * puedgergra * (1.0 / dx + 0.75 * met13uedger / dz)
        elseif ko + k == sizezz - nbz
            aru = 0.0
        else
            aru =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k + 1]) +
                guedger * puedgergra / dx
        end

        # ----------------- A(i+1,j,k-1) ------------------#

        if ko + k == k0
            ard = 0.0
        elseif ko + k == k0 + 1
            ard =
                -gedger * pedgergra * 0.25 * met13edger / dz +
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
                gdedger * pdedgergra * (1.0 / dx - 0.75 * met13dedger / dz)
        elseif ko + k == sizezz - nbz - 1
            ard =
                -gedger * pedgergra * 0.25 * met13edger / dz +
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
                gdedger * pdedgergra / dx +
                guedger * puedgergra * 0.25 * met13uedger / dz
        elseif ko + k == sizezz - nbz
            ard =
                -gedger * pedgergra * met13edger / dz +
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
                gdedger * pdedgergra / dx
        else
            ard =
                -gedger * pedgergra * 0.25 * met13edger / dz +
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i + 1, j, k] /
                (jac[i + 1, j, k] + jac[i + 1, j, k - 1]) +
                gdedger * pdedgergra / dx
        end

        # ----------------- A(i-1,j,k+1) ------------------#

        if ko + k == k0
            alu =
                gedgel * pedgelgra * met13edgel / dz -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                guedgel * puedgelgra / dx
        elseif ko + k == k0 + 1
            alu =
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                guedgel * puedgelgra / dx -
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        elseif ko + k == sizezz - nbz - 1
            alu =
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                guedgel * puedgelgra * (1.0 / dx - 0.75 * met13uedgel / dz)
        elseif ko + k == sizezz - nbz
            alu = 0.0
        else
            alu =
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgeu * pedgeugra * 0.5 * met13edgeu / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k + 1]) -
                guedgel * puedgelgra / dx
        end

        # ----------------- A(i-1,j,k-1) ------------------#

        if ko + k == k0
            ald = 0.0
        elseif ko + k == k0 + 1
            ald =
                -gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                gdedgel * pdedgelgra * (1.0 / dx + 0.75 * met13dedgel / dz)
        elseif ko + k == sizezz - nbz - 1
            ald =
                -gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                gdedgel * pdedgelgra / dx +
                guedgel * puedgelgra * 0.25 * met13uedgel / dz
        elseif ko + k == sizezz - nbz
            ald =
                -gedgel * pedgelgra * met13edgel / dz -
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                gdedgel * pdedgelgra / dx
        else
            ald =
                -gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedged * pedgedgra * 0.5 * met13edged / dx * jac[i - 1, j, k] /
                (jac[i - 1, j, k] + jac[i - 1, j, k - 1]) -
                gdedgel * pdedgelgra / dx
        end

        # ----------------- A(i,j+1,k+1) ------------------#

        if ko + k == k0
            afu =
                gedgef * pedgefgra * met23edgef / dz +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                guedgef * puedgefgra / dy
        elseif ko + k == k0 + 1
            afu =
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                guedgef * puedgefgra / dy -
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        elseif ko + k == sizezz - nbz - 1
            afu =
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                guedgef * puedgefgra * (1.0 / dy + 0.75 * met23uedgef / dz)
        elseif ko + k == sizezz - nbz
            afu = 0.0
        else
            afu =
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k + 1]) +
                guedgef * puedgefgra / dy
        end

        # ----------------- A(i,j+1,k-1) ------------------#

        if ko + k == k0
            afd = 0.0
        elseif ko + k == k0 + 1
            afd =
                -gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
                gdedgef * pdedgefgra * (1.0 / dy - 0.75 * met23dedgef / dz)
        elseif ko + k == sizezz - nbz - 1
            afd =
                -gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
                gdedgef * pdedgefgra / dy +
                guedgef * puedgefgra * 0.25 * met23uedgef / dz
        elseif ko + k == sizezz - nbz
            afd =
                -gedgef * pedgefgra * met23edgef / dz +
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
                gdedgef * pdedgefgra / dy
        else
            afd =
                -gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j + 1, k] /
                (jac[i, j + 1, k] + jac[i, j + 1, k - 1]) +
                gdedgef * pdedgefgra / dy
        end

        # ----------------- A(i,j-1,k+1) ------------------#

        if ko + k == k0
            abu =
                gedgeb * pedgebgra * met23edgeb / dz -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                guedgeb * puedgebgra / dy
        elseif ko + k == k0 + 1
            abu =
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                guedgeb * puedgebgra / dy -
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        elseif ko + k == sizezz - nbz - 1
            abu =
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                guedgeb * puedgebgra * (1.0 / dy - 0.75 * met23uedgeb / dz)
        elseif ko + k == sizezz - nbz
            abu = 0.0
        else
            abu =
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedgeu * pedgeugra * 0.5 * met23edgeu / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k + 1]) -
                guedgeb * puedgebgra / dy
        end

        # ----------------- A(i,j-1,k-1) ------------------#

        if ko + k == k0
            abd = 0.0
        elseif ko + k == k0 + 1
            abd =
                -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                gdedgeb * pdedgebgra * (1.0 / dy + 0.75 * met23dedgeb / dz)
        elseif ko + k == sizezz - nbz - 1
            abd =
                -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                gdedgeb * pdedgebgra / dy +
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k == sizezz - nbz
            abd =
                -gedgeb * pedgebgra * met23edgeb / dz -
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                gdedgeb * pdedgebgra / dy
        else
            abd =
                -gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gedged * pedgedgra * 0.5 * met23edged / dy * jac[i, j - 1, k] /
                (jac[i, j - 1, k] + jac[i, j - 1, k - 1]) -
                gdedgeb * pdedgebgra / dy
        end

        # ------------------ A(i,j,k+2) -------------------#

        if ko + k == k0
            auu =
                -gedger * pedgergra * 0.25 * met13edger / dz -
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gedgef * pedgefgra * 0.25 * met23edgef / dz -
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
                guedger * puedgergra * 0.25 * met13uedger / dz +
                guedgel * puedgelgra * 0.25 * met13uedgel / dz +
                guedgef * puedgefgra * 0.25 * met23uedgef / dz +
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k >= sizezz - nbz - 1
            auu = 0.0
        else
            auu =
                guedger * puedgergra * 0.25 * met13uedger / dz +
                guedgel * puedgelgra * 0.25 * met13uedgel / dz +
                guedgef * puedgefgra * 0.25 * met23uedgef / dz +
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        end

        # ------------------ A(i,j,k-2) -------------------#

        if ko + k <= k0 + 1
            add = 0.0
        elseif ko + k == sizezz - nbz
            add =
                gedger * pedgergra * 0.25 * met13edger / dz +
                gedgel * pedgelgra * 0.25 * met13edgel / dz +
                gedgef * pedgefgra * 0.25 * met23edgef / dz +
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gdedger * pdedgergra * 0.25 * met13dedger / dz -
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        else
            add =
                -gdedger * pdedgergra * 0.25 * met13dedger / dz -
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz -
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz -
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        end

        # ----------------- A(i+1,j,k+2) ------------------#

        if ko + k == k0
            aruu =
                -gedger * pedgergra * 0.25 * met13edger / dz +
                guedger * puedgergra * 0.25 * met13uedger / dz
        elseif ko + k >= sizezz - nbz - 1
            aruu = 0.0
        else
            aruu = guedger * puedgergra * 0.25 * met13uedger / dz
        end

        # ----------------- A(i+1,j,k-2) ------------------#

        if ko + k <= k0 + 1
            ardd = 0.0
        elseif ko + k == sizezz - nbz
            ardd =
                gedger * pedgergra * 0.25 * met13edger / dz -
                gdedger * pdedgergra * 0.25 * met13dedger / dz
        else
            ardd = -gdedger * pdedgergra * 0.25 * met13dedger / dz
        end

        # ----------------- A(i-1,j,k+2) ------------------#

        if ko + k == k0
            aluu =
                -gedgel * pedgelgra * 0.25 * met13edgel / dz +
                guedgel * puedgelgra * 0.25 * met13uedgel / dz
        elseif ko + k >= sizezz - nbz - 1
            aluu = 0.0
        else
            aluu = guedgel * puedgelgra * 0.25 * met13uedgel / dz
        end

        # ----------------- A(i-1,j,k-2) ------------------#

        if ko + k <= k0 + 1
            aldd = 0.0
        elseif ko + k == sizezz - nbz
            aldd =
                gedgel * pedgelgra * 0.25 * met13edgel / dz -
                gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        else
            aldd = -gdedgel * pdedgelgra * 0.25 * met13dedgel / dz
        end

        # ----------------- A(i,j+1,k+2) ------------------#

        if ko + k == k0
            afuu =
                -gedgef * pedgefgra * 0.25 * met23edgef / dz +
                guedgef * puedgefgra * 0.25 * met23uedgef / dz
        elseif ko + k >= sizezz - nbz - 1
            afuu = 0.0
        else
            afuu = guedgef * puedgefgra * 0.25 * met23uedgef / dz
        end

        # ----------------- A(i,j+1,k-2) ------------------#

        if ko + k <= k0 + 1
            afdd = 0.0
        elseif ko + k == sizezz - nbz
            afdd =
                gedgef * pedgefgra * 0.25 * met23edgef / dz -
                gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        else
            afdd = -gdedgef * pdedgefgra * 0.25 * met23dedgef / dz
        end

        # ----------------- A(i,j-1,k+2) ------------------#

        if ko + k == k0
            abuu =
                -gedgeb * pedgebgra * 0.25 * met23edgeb / dz +
                guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        elseif ko + k >= sizezz - nbz - 1
            abuu = 0.0
        else
            abuu = guedgeb * puedgebgra * 0.25 * met23uedgeb / dz
        end

        # ----------------- A(i,j-1,k-2) ------------------#

        if ko + k <= k0 + 1
            abdd = 0.0
        elseif ko + k == sizezz - nbz
            abdd =
                gedgeb * pedgebgra * 0.25 * met23edgeb / dz -
                gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        else
            abdd = -gdedgeb * pdedgebgra * 0.25 * met23dedgeb / dz
        end

        # Scale the tensor elements.
        ac = ac / (fcscal^2.0)
        ar = ar / fcscal / fcscal_r
        al = al / fcscal / fcscal_l
        af = af / fcscal / fcscal_f
        ab = ab / fcscal / fcscal_b
        au = au / fcscal / fcscal_u
        ad = ad / fcscal / fcscal_d
        aru = aru / fcscal / fcscal_ru
        ard = ard / fcscal / fcscal_rd
        alu = alu / fcscal / fcscal_lu
        ald = ald / fcscal / fcscal_ld
        afu = afu / fcscal / fcscal_fu
        afd = afd / fcscal / fcscal_fd
        abu = abu / fcscal / fcscal_bu
        abd = abd / fcscal / fcscal_bd
        auu = auu / fcscal / fcscal_uu
        add = add / fcscal / fcscal_dd
        aruu = aruu / fcscal / fcscal_ruu
        ardd = ardd / fcscal / fcscal_rdd
        aluu = aluu / fcscal / fcscal_luu
        aldd = aldd / fcscal / fcscal_ldd
        afuu = afuu / fcscal / fcscal_fuu
        afdd = afdd / fcscal / fcscal_fdd
        abuu = abuu / fcscal / fcscal_buu
        abdd = abdd / fcscal / fcscal_bdd

        # Determine indices for the operator.
        ia = i - i0 + 1
        ja = j - j0 + 1
        ka = k - k0 + 1

        # Set matrix elements for bicgstab.
        ac_b[ia, ja, ka] = ac
        ar_b[ia, ja, ka] = ar
        al_b[ia, ja, ka] = al
        af_b[ia, ja, ka] = af
        ab_b[ia, ja, ka] = ab
        au_b[ia, ja, ka] = au
        ad_b[ia, ja, ka] = ad
        aru_b[ia, ja, ka] = aru
        ard_b[ia, ja, ka] = ard
        alu_b[ia, ja, ka] = alu
        ald_b[ia, ja, ka] = ald
        afu_b[ia, ja, ka] = afu
        afd_b[ia, ja, ka] = afd
        abu_b[ia, ja, ka] = abu
        abd_b[ia, ja, ka] = abd
        auu_b[ia, ja, ka] = auu
        add_b[ia, ja, ka] = add
        aruu_b[ia, ja, ka] = aruu
        ardd_b[ia, ja, ka] = ardd
        aluu_b[ia, ja, ka] = aluu
        aldd_b[ia, ja, ka] = aldd
        afuu_b[ia, ja, ka] = afuu
        afdd_b[ia, ja, ka] = afdd
        abuu_b[ia, ja, ka] = abuu
        abdd_b[ia, ja, ka] = abdd
    end

    return
end
