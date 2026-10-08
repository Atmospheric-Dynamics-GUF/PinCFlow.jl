"""
    compute_spectral_cell_indices(
        state::State,
        kpr::AbstractFloat,
        mrr::AbstractFloat,
        dkpr::AbstractFloat,
        dmr::AbstractFloat,
    )::Union{Nothing, NTuple{4, <:Integer}}

From the given spectral ray-volume position and extent, determine the
indices of the spectral grid cells occupied by the ray volume.

For `x_size > 1`, a finite horizontal spectral extent `dkpr > 0`
may overlap multiple `kp` cells. If `dkpr == 0`, only the `kp` cell
containing the ray centre is returned.

For `x_size == 1`, the horizontal wavenumber is a discrete mode and
therefore belongs to exactly one `kp` index.

For `dmr > 0`, the vertical ray-volume extent may overlap multiple
`m` cells. If `dmr == 0`, only the `m` cell containing the ray centre
is returned.

If the ray volume lies completely outside the represented spectral
domain, `nothing` is returned.

For centre-based assignment, a point lying exactly on an internal
cell edge is assigned to the cell on the right. The uppermost domain
edge belongs to the final cell.
"""
function compute_spectral_cell_indices end

function compute_spectral_cell_indices(
    state::State,
    kpr::AbstractFloat,
    mrr::AbstractFloat,
    dkpr::AbstractFloat,
    dmr::AbstractFloat,
)::Union{Nothing, NTuple{4, <:Integer}}

    (; x_size) = state.namelists.domain
    (; kp, m, kpc, mc, kpl, ml) = state.spec_tend.spec_grid

    mrp = abs(mrr)
    m_l = mrp - dmr / 2
    m_u = mrp + dmr / 2

    kpmin = kpmax = mmin = mmax = 0

    # ------------------------------------------------------------------
    # Determine the positive-m part of the grid used for the bounds
    # and index search.
    # ------------------------------------------------------------------

    if m[1] < 0
        m0 = ml ÷ 2 + 1
        m1 = ml

        # First edge of the positive-m spectral domain.
        m_lo = mc[m0 + 1]
    else
        m0 = 1
        m1 = ml
        m_lo = mc[1]
    end

    m_hi = mc[end]

    # ------------------------------------------------------------------
    # Horizontal spectral indices.
    # ------------------------------------------------------------------

    if x_size == 1

        # k is a discrete Fourier mode:
        #
        #     kp[n] = n * 2π / Lx
        #
        # Therefore the ray belongs to exactly one kp index.
        kpi = compute_discrete_k_index(kp, kpr)

        if isnothing(kpi)
            return nothing
        end

        kpmin = kpi
        kpmax = kpi

    elseif iszero(dkpr)

        # --------------------------------------------------------------
        # Centre-based horizontal spectral assignment.
        #
        # The numerical ray-volume width does not participate in the
        # cell lookup. The ray centre belongs to exactly one kp cell.
        # --------------------------------------------------------------

        kp_lo = kpc[1]
        kp_hi = kpc[end]

        if kpr < kp_lo || kpr > kp_hi
            return nothing
        end

        for ii in 1:kpl
            if (
                kpc[ii] <= kpr < kpc[ii + 1] ||
                (ii == kpl && kpr == kpc[end])
            )
                kpmin = ii
                kpmax = ii
                break
            end
        end

    else

        # --------------------------------------------------------------
        # Existing finite-width horizontal spectral algorithm.
        # --------------------------------------------------------------

        kp_l = kpr - dkpr / 2
        kp_u = kpr + dkpr / 2

        kp_lo = kpc[1]
        kp_hi = kpc[end]

        # Ray volume lies completely outside the represented
        # horizontal spectral domain.
        if kp_u < kp_lo || kp_l > kp_hi
            return nothing
        end

        if kpl > 1
            for ii in eachindex(kp)
                if kpc[ii] <= kp_l
                    kpmin = ii
                end
                if kpc[ii + 1] < kp_u
                    kpmax = ii + 1
                end
            end
        end

        kpmin = clamp(kpmin, 1, kpl)
        kpmax = clamp(kpmax, 1, kpl)

    end

    # ------------------------------------------------------------------
    # Vertical spectral bounds.
    #
    # This remains applicable for both finite-width and centre-based
    # cases because for dmr == 0:
    #
    #     m_l = m_u = |mr|.
    # ------------------------------------------------------------------

    # Ray volume lies completely outside the represented vertical
    # spectral domain.
    if m_u < m_lo || m_l > m_hi
        return nothing
    end

    # ------------------------------------------------------------------
    # Vertical spectral indices.
    # ------------------------------------------------------------------

    if iszero(dmr)

        # --------------------------------------------------------------
        # Centre-based vertical spectral assignment.
        #
        # Work first with |m| on the positive-m half of the grid.
        # The result is then reflected back to the negative-m half
        # when mrr < 0.
        # --------------------------------------------------------------

        if m[1] < 0

            # Signed-m Triad2D grid.
            #
            # Positive-m cell jj is bounded by:
            #
            #     mc[jj + 1] <= |mr| < mc[jj + 2].
            #
            for jj in m0:m1
                if (
                    mc[jj + 1] <= mrp < mc[jj + 2] ||
                    (jj == m1 && mrp == mc[jj + 2])
                )
                    mmin = jj
                    mmax = jj
                    break
                end
            end

        else

            # Positive-m-only grid.
            #
            # Cell jj is bounded by:
            #
            #     mc[jj] <= mr < mc[jj + 1].
            #
            for jj in m0:m1
                if (
                    mc[jj] <= mrp < mc[jj + 1] ||
                    (jj == m1 && mrp == mc[jj + 1])
                )
                    mmin = jj
                    mmax = jj
                    break
                end
            end
        end

        # Convert positive-|m| index back to the negative-m half.
        if mrr < 0
            mmin = ml - mmin + 1
            mmax = mmin
        end

    else

        # --------------------------------------------------------------
        # Existing finite-width vertical spectral algorithm.
        # --------------------------------------------------------------

        if ml > 1
            if m[1] < 0
                # Signed-m Triad2D grid. Positive-m cell j is bounded by
                # mc[j + 1] and mc[j + 2].
                for jj in m0:m1
                    if mc[jj + 1] <= m_l
                        mmin = jj
                    end
                    if mc[jj + 2] < m_u
                        mmax = jj + 1
                    end
                end
            else
                # Positive-m-only grid.
                for jj in m0:m1
                    if mc[jj] <= m_l
                        mmin = jj
                    end
                    if mc[jj + 1] < m_u
                        mmax = jj + 1
                    end
                end
            end
        end

        mmin = clamp(mmin, m0, m1)
        mmax = clamp(mmax, m0, m1)

        # Convert positive-|m| indices back to the negative-m half.
        if mrr < 0
            mmax, mmin = ml - mmin + 1, ml - mmax + 1
        end
    end

    return (kpmin, kpmax, mmin, mmax)
end