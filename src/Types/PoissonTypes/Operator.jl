"""
    Operator{A <: AbstractArray{<:AbstractFloat, 3}}

Workspace array for matrix-vector operations in the Poisson solver.

# Fields

  - `s::A`: Auxiliary field for boundary communication and operator application

# Usage

Provides temporary storage for [`PinCFlow.PoissonSolver.apply_operator!`](@ref)
to handle boundary conditions and halo exchanges during matrix-vector products.
"""
struct Operator{A <: AbstractArray{<:AbstractFloat, 3}}
    s::A
end

"""
    Operator(domain::Domain)

Initialize operator workspace array sized according to extended domain.

# Arguments

  - `domain::Domain`: Domain specification with extended dimensions

# Returns

  - `Operator`: Container with zero-initialized workspace array
"""
function Operator(domain::Domain)

    # Get all necessary fields.
    (; nxx, nyy, nzz) = domain

    # Return an Operator instance.
    return Operator(zeros(nxx, nyy, nzz))
end
