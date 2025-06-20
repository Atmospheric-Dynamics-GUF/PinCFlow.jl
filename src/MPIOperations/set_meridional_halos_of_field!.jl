"""
    set_meridional_halos_of_field!(field::AbstractMatrix{<:AbstractFloat}, namelists::Namelists, domain::Domain)

Exchange meridional (y-direction) halo regions for 2D field arrays.

Performs bidirectional MPI communication between backward and forward neighbor
processes to maintain field continuity across y-direction domain boundaries.

# Arguments

  - `field::AbstractMatrix{<:AbstractFloat}`: 2D field array for halo exchange
  - `namelists::Namelists`: Configuration containing `nby` halo layer count
  - `domain::Domain`: MPI decomposition with neighbor process IDs

# Communication Pattern

  - **Send forward**: `field[:, (j1-nby+1):j1]` → forward neighbor's `(j0-nby):(j0-1)`
  - **Send backward**: `field[:, j0:(j0+nby-1)]` → backward neighbor's `(j1+1):(j1+nby)`
  - Uses `MPI.Sendrecv!` for deadlock-free bidirectional exchange
"""
function set_meridional_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; nby) = namelists.domain
    (; comm, j0, j1, backward, forward) = domain

    @views MPI.Sendrecv!(
        field[:, (j1 - nby + 1):j1],
        field[:, (j0 - nby):(j0 - 1)],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[:, j0:(j0 + nby - 1)],
        field[:, (j1 + 1):(j1 + nby)],
        comm;
        dest = backward,
        source = forward,
    )

    return
end

"""
    set_meridional_halos_of_field!(field::AbstractArray{<:Real, 3}, namelists::Namelists, domain::Domain; layers::NTuple{3, <:Integer} = (-1, -1, -1))

Exchange meridional (y-direction) halo regions for 3D field arrays.

# Arguments

  - `field::AbstractArray{<:Real, 3}`: 3D field array for halo exchange
  - `layers::NTuple{3, <:Integer}`: Custom halo sizes (nbx, nby, nbz). Use -1 for defaults

# Extended Region Coverage

  - **X-direction**: `(i0-nbx):(i1+nbx)` (includes zonal halos)
  - **Y-direction**: `nby` layers exchanged at each boundary
  - **Z-direction**: `(k0-nbz):(k1+nbz)` (includes vertical halos)
"""
function set_meridional_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, backward, forward) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[i, (j1 - nby + 1):j1, k],
        field[i, (j0 - nby):(j0 - 1), k],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[i, j0:(j0 + nby - 1), k],
        field[i, (j1 + 1):(j1 + nby), k],
        comm;
        dest = backward,
        source = forward,
    )

    return
end

"""
    set_meridional_halos_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists::Namelists, domain::Domain; layers::NTuple{3, <:Integer} = (-1, -1, -1))

Exchange meridional (y-direction) halo regions for 5D field arrays.

Handles multi-component fields such as metric tensors while preserving
tensor components in dimensions 4 and 5.

# Arguments

  - `field::AbstractArray{<:AbstractFloat, 5}`: 5D field array (e.g., metric tensors)
  - `layers::NTuple{3, <:Integer}`: Custom spatial halo sizes

# Tensor Handling

  - Exchanges spatial dimensions 1-3 between meridional neighbors
  - Preserves tensor components in dimensions 4-5
"""
function set_meridional_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, backward, forward) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[i, (j1 - nby + 1):j1, k, :, :],
        field[i, (j0 - nby):(j0 - 1), k, :, :],
        comm;
        dest = forward,
        source = backward,
    )

    @views MPI.Sendrecv!(
        field[i, j0:(j0 + nby - 1), k, :, :],
        field[i, (j1 + 1):(j1 + nby), k, :, :],
        comm;
        dest = backward,
        source = forward,
    )

    return
end
