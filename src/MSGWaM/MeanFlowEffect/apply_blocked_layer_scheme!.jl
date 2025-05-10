function apply_blocked_layer_scheme!(state::State)
    (; testcase) = state.namelists.setting
    apply_blocked_layer_scheme!(state, testcase)
    return
end

function apply_blocked_layer_scheme!(
    state::State,
    testcase::AbstractWKBTestCase,
)
    return
end

function apply_blocked_layer_scheme!(state::State, testcase::WKBMountainWave)
    (; blocking, drag_coefficient) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain
    (; dz, jac, ztildetfc, topography_spectrum, k_spectrum, l_spectrum) =
        state.grid
    (; rhostrattfc) = state.atmosphere
    (; rho, u, v) = state.variables.predictands
    (; zb) = state.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    if !blocking
        return
    end

    # Initialize arrays for blocked-flow drag computation.
    (kavg, uperp, drag) = (zeros(2) for i in 1:3)

    # Adjust the drag to account for blocking.
    for kz in k0:k1, jy in j0:j1, ix in i0:i1
        fraction =
            (
                min(zb[ix, jy], ztildetfc[ix, jy, kz]) -
                ztildetfc[ix, jy, kz - 1]
            ) / jac[ix, jy, kz] / dz
        @views if fraction <= 0
            continue
        else
            kavg[1] =
                sum(
                    abs.(topography_spectrum[:, ix, jy]) .*
                    k_spectrum[:, ix, jy],
                ) / sum(abs.(topography_spectrum[:, ix, jy]))
            kavg[2] =
                sum(
                    abs.(topography_spectrum[:, ix, jy]) .*
                    l_spectrum[:, ix, jy],
                ) / sum(abs.(topography_spectrum[:, ix, jy]))

            uperp .=
                (
                    (u[ix, jy, kz] .+ u[ix - 1, jy, kz]) .* kavg[1] .+
                    (v[ix, jy, kz] .+ v[ix, jy - 1, kz]) .* kavg[2]
                ) ./ 2 .* kavg ./ dot(kavg, kavg)
            drag .=
                -drag_coefficient .*
                (rho[ix, jy, kz] .+ rhostrattfc[ix, jy, kz]) .*
                sqrt(dot(kavg, kavg)) ./ (2 .* pi) .* sqrt(dot(uperp, uperp)) .*
                uperp
            dudt[ix, jy, kz] =
                fraction * drag[1] + (1 - fraction) * dudt[ix, jy, kz]
            dvdt[ix, jy, kz] =
                fraction * drag[2] + (1 - fraction) * dvdt[ix, jy, kz]
            dthetadt[ix, jy, kz] = (1 - fraction) * dthetadt[ix, jy, kz]
        end
    end
end
