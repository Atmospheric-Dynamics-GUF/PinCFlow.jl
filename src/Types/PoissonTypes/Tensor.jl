"""
    Tensor{A <: AbstractArray{<:AbstractFloat, 3}}

Storage for 27-point stencil matrix coefficients of the Poisson operator.

# Fields

  - `ac_b::A`: Center coefficient
  - `acv_b::A`: Center vertical coefficient
  - `ach_b::A`: Center horizontal coefficient
  - `al_b::A`: Left neighbor coefficient
  - `ar_b::A`: Right neighbor coefficient
  - `ab_b::A`: Back neighbor coefficient
  - `af_b::A`: Front neighbor coefficient
  - `ad_b::A`: Down neighbor coefficient
  - `au_b::A`: Up neighbor coefficient
  - `ald_b::A`: Left-down diagonal coefficient
  - `alu_b::A`: Left-up diagonal coefficient
  - `ard_b::A`: Right-down diagonal coefficient
  - `aru_b::A`: Right-up diagonal coefficient
  - `abd_b::A`: Back-down diagonal coefficient
  - `abu_b::A`: Back-up diagonal coefficient
  - `afd_b::A`: Front-down diagonal coefficient
  - `afu_b::A`: Front-up diagonal coefficient
  - `add_b::A`: Down-down coefficient
  - `auu_b::A`: Up-up coefficient
  - `aldd_b::A`: Left-down-down coefficient
  - `aluu_b::A`: Left-up-up coefficient
  - `ardd_b::A`: Right-down-down coefficient
  - `aruu_b::A`: Right-up-up coefficient
  - `abdd_b::A`: Back-down-down coefficient
  - `abuu_b::A`: Back-up-up coefficient
  - `afdd_b::A`: Front-down-down coefficient
  - `afuu_b::A`: Front-up-up coefficient

# Usage

Matrix coefficients are computed by [`PinCFlow.PoissonSolver.compute_operator!`](@ref)
and used by [`PinCFlow.PoissonSolver.apply_operator!`](@ref) for matrix-vector products.
"""
struct Tensor{A <: AbstractArray{<:AbstractFloat, 3}}
    ac_b::A
    acv_b::A
    ach_b::A
    al_b::A
    ar_b::A
    ab_b::A
    af_b::A
    ad_b::A
    au_b::A
    ald_b::A
    alu_b::A
    ard_b::A
    aru_b::A
    abd_b::A
    abu_b::A
    afd_b::A
    afu_b::A
    add_b::A
    auu_b::A
    aldd_b::A
    aluu_b::A
    ardd_b::A
    aruu_b::A
    abdd_b::A
    abuu_b::A
    afdd_b::A
    afuu_b::A
end

"""
    Tensor(domain::Domain)

Initialize tensor coefficient arrays sized according to local domain.

# Arguments

  - `domain::Domain`: Local domain dimensions

# Returns

  - `Tensor`: Container with zero-initialized coefficient arrays
"""
function Tensor(domain::Domain)
    # Get parameters.
    (; nx, ny, nz) = domain

    # Return a Tensor instance.
    return Tensor([zeros(nx, ny, nz) for i in 1:27]...)
end
