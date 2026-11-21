
# ============================================================
# Generic prescribed lower-boundary source interface
#
# Source properties and physical parameters are defined
# externally and passed as a callable object.
#
# No assumption is made about the spatial or temporal
# dependence of the source.
# ============================================================

"""
Set the physical time used when evaluating the boundary source.
"""
@inline function set_boundary_time!(state::State, time_seconds::Real)
    state.boundary_time[] = Float64(time_seconds)
    return nothing
end

"""
Evaluate the externally supplied boundary source.

The source must be callable as

    source(x, y, z, t)

where x, y, z and t are dimensional coordinates and time.

The source is free to ignore any of these arguments.
"""
@inline function evaluate_boundary_source(
    state::State,
    xdim::Real,
    ydim::Real,
    zdim::Real,
)
    source = state.boundary_source

    source === nothing &&
        error("No boundary source has been provided")

    return source(xdim, ydim, zdim, state.boundary_time[])
end
