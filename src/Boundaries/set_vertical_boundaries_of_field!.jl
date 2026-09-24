"""
```julia
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
```

Enforce vertical boundary conditions by dispatching to the vertical boundary condition appropriate method.

```julia 
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function,
    vertical_boundary_condition::Val{:Periodic};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
```

Enforce periodic vertical boundary conditions for a 3D array.

Halo exchange is used for multi-process domains (`npz > 1`).

```julia 
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function,
    vertical_boundary_condition::Val{:SolidWall};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
```

Enforce solid-wall vertical boundary conditions for a 3D array.

Halo exchange is used for multi-process domains (`npz > 1`). Use `mode = +` (`mode = -`) for line-reflected (point-reflected) ghost-cell values.

```julia 
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce vertical boundary conditions for a 5D array by dispatching to the appropriate vertical boundary condition method.

```julia
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    vertical_boundary_condition::Val{:Periodic};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Enforce periodic vertical boundary conditions for a 5D array. 

Halo exchange is used for multi-process domains (`npz > 1`).

```julia
set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
```

Exchange halo values of a 5D array if multiple processes are used in the vertical (`npz > 1`).

This method is applied to reconstruction arrays. Vertical boundary conditions are not enforced for these but for the fluxes determined from them.

# Arguments

  - `field`: Input array.

  - `namelists`: Namelists with all model parameters.

  - `domain`: Collection of domain-decomposition and MPI-communication parameters.

  - `mode`: Method used for setting the boundary-cell values.

  - `vertical_boundary_condition`: Vertical boundary condition.

# Keywords

  - `layers`: The number of boundary layers in each dimension. Use `-1` for the default values from `namelists`.

  - `staggered`: A switch for whether or not the field is on the staggered vertical grid.

# See also

  - [`PinCFlow.MPIOperations.set_vertical_halos_of_field!`](@ref)
"""
function set_vertical_boundaries_of_field! end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
    (; vertical_boundary_condition) = namelists.domain

    @dispatch_vertical_boundary_condition set_vertical_boundaries_of_field!(
        field,
        namelists,
        domain,
        mode,
        Val(vertical_boundary_condition);
        layers,
        staggered,
    )

    return
end

@ivy function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function,
    vertical_boundary_condition::Val{:Periodic};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
    (; z_size) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if z_size > 1
        set_vertical_halos_of_field!(field, namelists, domain; layers)
    else
        ii = (i0 - nbx):(i1 + nbx)
        jj = (j0 - nby):(j1 + nby)

        for k in 1:nbz
            field[ii, jj, k0 - k] .= field[ii, jj, k1 - k + 1]
            field[ii, jj, k1 + k] .= field[ii, jj, k0 + k - 1]
        end
    end

    return
end

@ivy function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:Real, 3},
    namelists::Namelists,
    domain::Domain,
    mode::Function,
    vertical_boundary_condition::Val{:SolidWall};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
    staggered = false,
)
    (; z_size) = namelists.domain
    (; nz, ko, i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if z_size > 1
        set_vertical_halos_of_field!(field, namelists, domain; layers)
    end

    ii = (i0 - nbx):(i1 + nbx)
    jj = (j0 - nby):(j1 + nby)

    if ko == 0
        if staggered
            field[ii, jj, k0 - 1] .= 0.0
            for k in 1:nbz
                field[ii, jj, k0 - k] .= mode.(field[ii, jj, k0 + k - 2])
            end
        else
            for k in 1:nbz
                field[ii, jj, k0 - k] .= mode.(field[ii, jj, k0 + k - 1])
            end
        end
    end

    if ko + nz == z_size
        if staggered
            field[ii, jj, k1] .= 0.0
            for k in 1:nbz
                field[ii, jj, k1 + k] .= mode.(field[ii, jj, k1 - k])
            end
        else
            for k in 1:nbz
                field[ii, jj, k1 + k] .= mode.(field[ii, jj, k1 - k + 1])
            end
        end
    end

    return
end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain;
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; vertical_boundary_condition) = namelists.domain

    @dispatch_vertical_boundary_condition set_vertical_boundaries_of_field!(
        field,
        namelists,
        domain,
        Val(vertical_boundary_condition);
        layers,
    )
    return
end

@ivy function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    vertical_boundary_condition::Val{:Periodic};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; z_size) = namelists.domain
    (; i0, i1, j0, j1, k0, k1) = domain

    nbx = layers[1] == -1 ? namelists.domain.nbx : layers[1]
    nby = layers[2] == -1 ? namelists.domain.nby : layers[2]
    nbz = layers[3] == -1 ? namelists.domain.nbz : layers[3]

    if z_size > 1
        set_vertical_halos_of_field!(field, namelists, domain; layers)
    else
        ii = (i0 - nbx):(i1 + nbx)
        jj = (j0 - nby):(j1 + nby)

        for k in 1:nbz
            field[ii, jj, k0 - k, :, :] .= field[ii, jj, k1 - k + 1, :, :]
            field[ii, jj, k1 + k, :, :] .= field[ii, jj, k0 + k - 1, :, :]
        end
    end

    return
end

function set_vertical_boundaries_of_field!(
    field::AbstractArray{<:AbstractFloat, 5},
    namelists::Namelists,
    domain::Domain,
    vertical_boundary_condition::Val{:SolidWall};
    layers::NTuple{3, <:Integer} = (-1, -1, -1),
)
    (; z_size) = namelists.domain

    if z_size > 1
        set_vertical_halos_of_field!(field, namelists, domain; layers)
    end

    return
end
