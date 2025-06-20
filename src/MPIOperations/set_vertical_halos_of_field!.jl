"""
    set_vertical_halos_of_field!(field::AbstractArray{<:Real, 3}, namelists::Namelists, domain::Domain, zboundaries::SolidWallBoundaries; layers::NTuple{3, <:Integer} = (-1, -1, -1))

Exchange vertical (z-direction) halo regions for 3D field arrays with solid wall boundaries.

Performs MPI communication to update ghost/halo cells at the top and bottom
boundaries of the local domain. Handles special cases for processes at the
global domain boundaries with solid wall boundary conditions.

# Arguments

  - `field::AbstractArray{<:Real, 3}`: 3D field array to exchange halos for
  - `namelists::Namelists`: Configuration containing boundary layer sizes
  - `domain::Domain`: MPI domain decomposition information
  - `zboundaries::SolidWallBoundaries`: Solid wall boundary condition type
  - `layers::NTuple{3, <:Integer}`: Custom halo sizes (nbx, nby, nbz). Use -1 for default from namelists

# Boundary Handling

  - **Bottom boundary process** (`ko == 0`): Only exchanges with upper neighbor
  - **Top boundary process** (`ko + nzz == sizezz`): Only exchanges with lower neighbor
  - **Interior processes**: Exchange with both upper and lower neighbors

# Communication Pattern

  - Sends top `nbz` layers to upper neighbor, receives from lower neighbor
  - Sends bottom `nbz` layers to lower neighbor, receives from upper neighbor
  - Includes extended regions in x and y directions for complete halo coverage

# Solid Wall Treatment

At solid wall boundaries, the halo exchange respects the no-penetration condition
by appropriately handling the boundary layer communications.
"""
function set_vertical_halos_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1],
            field[i, j, (k1 + 1):(k1 + nbz)],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1)],
            field[i, j, (k0 - nbz):(k0 - 1)],
            comm;
            dest = down,
            source = down,
        )
    else
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1],
            field[i, j, (k0 - nbz):(k0 - 1)],
            comm;
            dest = up,
            source = down,
        )

        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1)],
            field[i, j, (k1 + 1):(k1 + nbz)],
            comm;
            dest = down,
            source = up,
        )
    end

    return
end

"""
    set_vertical_halos_of_field!(field::AbstractArray{<:AbstractFloat, 5}, namelists::Namelists, domain::Domain, zboundaries::SolidWallBoundaries; layers::NTuple{3, <:Integer} = (-1, -1, -1))

Exchange vertical (z-direction) halo regions for 5D field arrays with solid wall boundaries.

Performs MPI communication to update ghost/halo cells for 5D arrays such as
metric tensors, handling solid wall boundary conditions at the domain top/bottom.

# Arguments

  - `field::AbstractArray{<:AbstractFloat, 5}`: 5D field array to exchange halos for
  - `namelists::Namelists`: Configuration containing boundary layer sizes
  - `domain::Domain`: MPI domain decomposition information
  - `zboundaries::SolidWallBoundaries`: Solid wall boundary condition type
  - `layers::NTuple{3, <:Integer}`: Custom halo sizes (nbx, nby, nbz). Use -1 for default from namelists

# Boundary Handling

  - **Bottom boundary process** (`ko == 0`): Exchanges only with upper neighbor
  - **Top boundary process** (`ko + nzz == sizezz`): Exchanges only with lower neighbor
  - **Interior processes**: Full bidirectional exchange with upper and lower neighbors

# Communication Pattern

  - Exchanges first 3 dimensions (spatial) between vertical neighbors
  - Preserves all data in dimensions 4 and 5 (e.g., tensor components)
  - Uses extended x and y regions to maintain corner/edge consistency

# Use Cases

Typically used for:

  - Metric tensor fields requiring vertical derivatives
  - 5D arrays in terrain-following coordinate systems
  - Multi-component fields near solid boundaries
"""
function set_vertical_halos_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    zboundaries::SolidWallBoundaries;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; nbz) = namelists.domain
    (; comm, sizezz, nzz, ko, i0, i1, j0, j1, k0, k1, down, up) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    i = (i0 - nbx):(i1 + nbx)
    j = (j0 - nby):(j1 + nby)

    if ko == 0
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1, :, :],
            field[i, j, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = up,
            source = up,
        )
    elseif ko + nzz == sizezz
        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1), :, :],
            field[i, j, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = down,
            source = down,
        )
    else
        @views MPI.Sendrecv!(
            field[i, j, (k1 - nbz + 1):k1, :, :],
            field[i, j, (k0 - nbz):(k0 - 1), :, :],
            comm;
            dest = up,
            source = down,
        )

        @views MPI.Sendrecv!(
            field[i, j, k0:(k0 + nbz - 1), :, :],
            field[i, j, (k1 + 1):(k1 + nbz), :, :],
            comm;
            dest = down,
            source = up,
        )
    end

    return
end
