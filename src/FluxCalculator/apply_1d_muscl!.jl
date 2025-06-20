"""
    apply_1d_muscl!(phi, phitilde, phisize)

Applies the Monotonic Upstream-centered Scheme for Conservation Laws (MUSCL) reconstruction
in one dimension.

# Arguments

  - `phi::AbstractVector{<:AbstractFloat}`: Input scalar values at cell centers.
  - `phitilde::AbstractMatrix{<:AbstractFloat}`: Output reconstructed values at cell interfaces.
    Column 1 contains (left,down,back) interface values, column 2 contains (right,up,front) interface values.
  - `phisize::Integer`: Size of the input array `phi`, including ghost cells.

# TODO: this assumes that phi has only one ghost cell on each side?

# Notes

  - The function modifies `phitilde` in-place.
  - Boundary cells (1 and phisize) are not reconstructed.
"""
function apply_1d_muscl!(
    phi::AbstractVector{<:AbstractFloat},
    phitilde::AbstractMatrix{<:AbstractFloat},
    phisize::Integer,
)

    # Initialize phitilde.
    phitilde .= 1000.0

    # Reconstruct.
    for i in 2:(phisize - 1)
        deltal = phi[i] - phi[i - 1]
        deltar = phi[i + 1] - phi[i]

        if deltar == 0.0
            phitilde[i, 2] = phi[i]
            phitilde[i, 1] = phi[i]
        else
            if deltal == 0.0
                theta = deltal / deltar
                s = (2.0 + theta) / 3.0
                sigmal = max(0.0, min(2 * theta, s, 2.0))

                phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
                phitilde[i, 1] = phi[i]
            else
                theta = deltal / deltar

                s = (2.0 + theta) / 3.0
                sigmal = max(0.0, min(2 * theta, s, 2.0))

                s = (2.0 + 1.0 / theta) / 3.0
                sigmar = max(0.0, min(2 / theta, s, 2.0))

                phitilde[i, 2] = phi[i] + 0.5 * sigmal * deltar
                phitilde[i, 1] = phi[i] - 0.5 * sigmar * deltal
            end
        end
    end

    # Return.
    return
end
