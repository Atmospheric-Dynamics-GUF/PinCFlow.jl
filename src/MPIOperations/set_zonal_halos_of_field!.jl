"""
```julia
set_zonal_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
```

Exchange zonal (`x`-direction) halo regions for 2D field arrays.

Performs bidirectional MPI communication between left and right neighbor
processes to maintain field continuity across domain boundaries in `x`-direction.

# Arguments

  - `field::AbstractMatrix{<:AbstractFloat}`: 2D field array to exchange halos for
  - `namelists::Namelists`: Configuration containing boundary layer sizes
  - `domain::Domain`: MPI domain decomposition information

# Communication Pattern

  - Sends rightmost `nbx` columns to right neighbor, receives into left halos
  - Sends leftmost `nbx` columns to left neighbor, receives into right halos
  - Uses `MPI.Sendrecv!` for bidirectional communication
"""
function set_zonal_halos_of_field!(
    field::AbstractMatrix{<:AbstractFloat},
    namelists::Namelists,
    domain::Domain,
)
    (; nbx) = namelists.domain
    (; comm, i0, i1, left, right) = domain

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, :],
        field[(i0 - nbx):(i0 - 1), :],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), :],
        field[(i1 + 1):(i1 + nbx), :],
        comm;
        dest = left,
        source = right,
    )

    return
end

"""
```julia
set_zonal_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange zonal (x-direction) halo regions for 3D field arrays.

Performs MPI communication to update ghost/halo cells at the left and right
boundaries of the local domain in the x-direction for 3D arrays.

# Arguments

  - `field::AbstractArray{<:Real, 3}`: 3D field array to exchange halos for
  - `namelists::Namelists`: Configuration containing boundary layer sizes
  - `domain::Domain`: MPI domain decomposition information
  - `layers::NTuple{3, <:Integer}`: Custom halo sizes (nbx, nby, nbz). Use -1 for default from namelists

# Communication Pattern

  - Exchanges data in x-direction between left and right neighbor processes
  - Includes extended regions in y and z directions to maintain corner/edge halos
  - Uses `MPI.Sendrecv!` for simultaneous send/receive operations

# Details

The exchange includes:

  - x-direction: `nbx` layers on each side
  - y-direction: extended by `nby` layers on each side
  - z-direction: extended by `nbz` layers on each side
"""
function set_zonal_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, j, k],
        field[(i0 - nbx):(i0 - 1), j, k],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), j, k],
        field[(i1 + 1):(i1 + nbx), j, k],
        comm;
        dest = left,
        source = right,
    )

    return
end

"""
```julia
set_zonal_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange zonal (x-direction) halo regions for 5D field arrays.

Performs MPI communication to update ghost/halo cells for 5D arrays such as
metric tensors or multi-component fields.

# Arguments

  - `field::AbstractArray{<:AbstractFloat, 5}`: 5D field array to exchange halos for
  - `namelists::Namelists`: Configuration containing boundary layer sizes
  - `domain::Domain`: MPI domain decomposition information
  - `layers::NTuple{3, <:Integer}`: Custom halo sizes (nbx, nby, nbz). Use -1 for default from namelists

# Communication Pattern

  - Exchanges first 3 dimensions (spatial) between left and right neighbors
  - Preserves all data in dimensions 4 and 5 (e.g., tensor components)
  - Uses extended regions in y and z directions for complete halo coverage

# Use Cases

  - Metric tensor fields with shape (nx, ny, nz, 3, 3)
  - Multi-component vector/tensor fields
  - Fields requiring spatial derivatives near boundaries
"""
function set_zonal_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; comm, i0, i1, j0, j1, k0, k1, left, right) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    j = (j0 - nby):(j1 + nby)
    k = (k0 - nbz):(k1 + nbz)

    @views MPI.Sendrecv!(
        field[(i1 - nbx + 1):i1, j, k, :, :],
        field[(i0 - nbx):(i0 - 1), j, k, :, :],
        comm;
        dest = right,
        source = left,
    )

    @views MPI.Sendrecv!(
        field[i0:(i0 + nbx - 1), j, k, :, :],
        field[(i1 + 1):(i1 + nbx), j, k, :, :],
        comm;
        dest = left,
        source = right,
    )

    return
end
