"""
    copy_rays!(rays::Rays, source::NTuple{4, <:Integer}, target::NTuple{4, <:Integer})

Copy all properties of a ray volume from source to target location.

Performs a complete copy of all ray volume properties from one array location
to another, including physical positions, spectral properties, extents, and
wave action density.

# Arguments

  - `rays::Rays`: Ray volume data structure containing all ray arrays
  - `source::NTuple{4, <:Integer}`: Source indices (iray_src, ix_src, jy_src, kz_src)
  - `target::NTuple{4, <:Integer}`: Target indices (iray_tgt, ix_tgt, jy_tgt, kz_tgt)

# Copied Properties

Copies all fields of the `Rays` structure:

  - **Physical position**: `x`, `y`, `z`
  - **Spectral position**: `k`, `l`, `m` (wavenumbers)
  - **Physical extents**: `dxray`, `dyray`, `dzray`
  - **Spectral extents**: `dkray`, `dlray`, `dmray`
  - **Wave action density**: `dens`

# Use Cases

## Ray Shifting

When rays move between grid cells:

```julia
copy_rays!(rays, (iray, ix_old, jy, kz), (jray, ix_new, jy, kz))
```

## Ray Compaction

Removing gaps in ray arrays:

```julia
copy_rays!(rays, (iray_sparse, ix, jy, kz), (iray_compact, ix, jy, kz))
```

## Ray Splitting

Duplicating rays before modification:

```julia
copy_rays!(rays, (iray, ix, jy, kz), (iray_new, ix, jy, kz))
```

## Boundary Conditions

Transferring rays to halo regions:

```julia
copy_rays!(rays, (iray, ix, jy, kz), (jray, ix_halo, jy, kz))
```

# Implementation

Uses reflection to iterate over all fields in the `Rays` structure,
ensuring complete copying even if new fields are added to the structure.

# Memory Considerations

  - Does not allocate new memory (in-place copy)
  - Source and target can be the same location (no-op)
  - Efficient field-by-field assignment

# Thread Safety

Safe for concurrent read access to source, but target location should
not be accessed concurrently during the copy operation.
"""
function copy_rays!(
    rays::Rays,
    source::NTuple{4, <:Integer},
    target::NTuple{4, <:Integer},
)
    (irs, ixs, jys, kzs) = source
    (irt, ixt, jyt, kzt) = target

    for field in fieldnames(Rays)
        getfield(rays, field)[irt, ixt, jyt, kzt] =
            getfield(rays, field)[irs, ixs, jys, kzs]
    end

    return
end
