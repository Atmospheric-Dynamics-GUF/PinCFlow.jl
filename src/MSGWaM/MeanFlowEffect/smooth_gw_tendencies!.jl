"""
    smooth_gw_tendencies!(state::State)

Apply spatial smoothing to gravity wave tendency fields.

Smooths the computed gravity wave momentum and heating tendencies to
reduce noise and improve numerical stability, using the specified
filter type and domain dimensionality.

# Arguments

  - `state::State`: Complete simulation state containing GW tendencies

# Domain-Dependent Smoothing

Applies smoothing only in active spatial dimensions:

  - **1D (Column)**: Z-direction only
  - **2D (XZ-plane)**: X and Z directions
  - **2D (YZ-plane)**: Y and Z directions
  - **3D (Full)**: X, Y, and Z directions

# Smoothed Fields

  - `dudt`: Zonal wind tendency
  - `dvdt`: Meridional wind tendency
  - `dthetadt`: Potential temperature tendency

# Filter Selection

Dispatches to appropriate filter based on `sm_filter` setting:

  - `Box`: Simple box filter (moving average)
  - `Shapiro`: Shapiro filter (selective smoothing)
"""
function smooth_gw_tendencies!(state::State)
    (; sizex, sizey) = state.namelists.domain
    (; lsmth_wkb, sm_filter) = state.namelists.wkb
    (; dudt, dvdt, dthetadt) = state.wkb.tendencies

    if !lsmth_wkb
        return
    end

    if sizex == sizey == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, Z())
        smooth_gw_tendencies!(dvdt, state, sm_filter, Z())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, Z())
    elseif sizex == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, YZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, YZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, YZ())
    elseif sizey == 1
        smooth_gw_tendencies!(dudt, state, sm_filter, XZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, XZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, XZ())
    else
        smooth_gw_tendencies!(dudt, state, sm_filter, XYZ())
        smooth_gw_tendencies!(dvdt, state, sm_filter, XYZ())
        smooth_gw_tendencies!(dthetadt, state, sm_filter, XYZ())
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Box, direction::XYZ)

Apply 3D box filter smoothing in all spatial directions.

Applies simple moving average smoothing using a cubic box filter
in all three spatial dimensions.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state containing grid and smoothing parameters
  - `sm_filter::Box`: Box filter type
  - `direction::XYZ`: 3D smoothing direction

# Algorithm

For each interior point, replaces value with average of
`(2*nsmth_wkb + 1)³` surrounding points.

# Boundary Requirements

Requires sufficient halo points: `nbx, nby, nbz ≥ nsmth_wkb`
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XYZ,
)
    (; nbx, nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end
    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    (i - nsmth_wkb):(i + nsmth_wkb),
                    (j - nsmth_wkb):(j + nsmth_wkb),
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^3
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Box, direction::XZ)

Apply 2D box filter smoothing in X and Z directions.

Applies moving average smoothing in the XZ-plane, leaving
Y-direction unsmoothed for 2D or column simulations.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Box`: Box filter type
  - `direction::XZ`: 2D smoothing in XZ-plane

# Algorithm

For each point, averages over `(2*nsmth_wkb + 1)²` points
in the XZ-plane at each Y-level.
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::XZ,
)
    (; nbx, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    (i - nsmth_wkb):(i + nsmth_wkb),
                    j,
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^2
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Box, direction::YZ)

Apply 2D box filter smoothing in Y and Z directions.

Similar to XZ smoothing but operates in the YZ-plane.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Box`: Box filter type
  - `direction::YZ`: 2D smoothing in YZ-plane
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::YZ,
)
    (; nby, nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end
    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(
                input[
                    i,
                    (j - nsmth_wkb):(j + nsmth_wkb),
                    (k - nsmth_wkb):(k + nsmth_wkb),
                ],
            ) / (2 * nsmth_wkb + 1)^2
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Box, direction::Z)

Apply 1D box filter smoothing in Z direction only.

Applies vertical smoothing for single-column or when only
vertical filtering is desired.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Box`: Box filter type
  - `direction::Z`: 1D smoothing in Z-direction

# Algorithm

For each horizontal point, averages over `2*nsmth_wkb + 1`
vertical levels.
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Box,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; i0, i1, j0, j1, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for k in k0:k1, j in j0:j1, i in i0:i1
        @views output[i, j, k] =
            sum(input[i, j, (k - nsmth_wkb):(k + nsmth_wkb)]) /
            (2 * nsmth_wkb + 1)
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::XYZ)

Apply 3D Shapiro filter smoothing in all spatial directions.

Applies Shapiro filter sequentially in X, Y, and Z directions.
Shapiro filters provide selective smoothing that preserves
sharp gradients better than box filters.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::XYZ`: 3D smoothing direction

# Sequential Application

 1. Apply Shapiro filter in X-direction
 2. Apply Shapiro filter in Y-direction
 3. Apply Shapiro filter in Z-direction
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XYZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, X())
    smooth_gw_tendencies!(output, state, sm_filter, Y())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::XZ)

Apply 2D Shapiro filter smoothing in X and Z directions.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::XZ`: 2D smoothing in XZ-plane
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::XZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, X())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::YZ)

Apply 2D Shapiro filter smoothing in Y and Z directions.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::YZ`: 2D smoothing in YZ-plane
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::YZ,
)
    smooth_gw_tendencies!(output, state, sm_filter, Y())
    smooth_gw_tendencies!(output, state, sm_filter, Z())
    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::Z)

Apply 1D Shapiro filter smoothing in Z direction.

Uses `apply_shapiro_filter!` function to apply selective smoothing
in the vertical direction for each horizontal column.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::Z`: 1D smoothing in Z-direction
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Z,
)
    (; nbz) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nxx, nyy, k0, k1) = state.domain

    if nbz < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbz < nsmth_wkb!")
    end

    input = copy(output)
    for j in 1:nyy, i in 1:nxx
        @views apply_shapiro_filter!(
            output[i, j, :],
            input[i, j, :],
            (k0, k1),
            Val(nsmth_wkb),
        )
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::Y)

Apply 1D Shapiro filter smoothing in Y direction.

Applies Shapiro filter along meridional direction for each
vertical column in the XZ-plane.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::Y`: 1D smoothing in Y-direction
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::Y,
)
    (; nby) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nxx, nzz, j0, j1) = state.domain

    if nby < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nby < nsmth_wkb!")
    end

    input = copy(output)
    for k in 1:nzz, i in 1:nxx
        @views apply_shapiro_filter!(
            output[i, :, k],
            input[i, :, k],
            (j0, j1),
            Val(nsmth_wkb),
        )
    end

    return
end

"""
    smooth_gw_tendencies!(output::AbstractArray{<:AbstractFloat, 3}, state::State, sm_filter::Shapiro, direction::X)

Apply 1D Shapiro filter smoothing in X direction.

Applies Shapiro filter along zonal direction for each
vertical column in the YZ-plane.

# Arguments

  - `output`: Field to smooth (modified in-place)
  - `state`: Simulation state
  - `sm_filter::Shapiro`: Shapiro filter type
  - `direction::X`: 1D smoothing in X-direction
"""
function smooth_gw_tendencies!(
    output::AbstractArray{<:AbstractFloat, 3},
    state::State,
    sm_filter::Shapiro,
    direction::X,
)
    (; nbx) = state.namelists.domain
    (; nsmth_wkb) = state.namelists.wkb
    (; nyy, nzz, i0, i1) = state.domain

    if nbx < nsmth_wkb
        error("Error in smooth_gw_tendencies!: nbx < nsmth_wkb!")
    end

    input = copy(output)
    for k in 1:nzz, j in 1:nyy
        @views apply_shapiro_filter!(
            output[:, j, k],
            input[:, j, k],
            (i0, i1),
            Val(nsmth_wkb),
        )
    end

    return
end
