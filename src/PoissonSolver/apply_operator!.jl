"""
```julia
apply_operator!(
    sin::AbstractArray{<:AbstractFloat, 3},
    ls::AbstractArray{<:AbstractFloat, 3},
    hortot::Total,
    state::State,
)
```

Apply the total linear operator to the solution array `sin`.

Before the operator is applied, the boundary/halo values of `sin` are set, using `set_zonal_boundaries_of_field!`, `set_meridional_boundaries_of_field!` and `set_vertical_halos_of_field!`. Note that in the vertical, only halo values need to be set (if the domain is parallelized in that direction), due to the solid-wall boundaries.

```julia
apply_operator!(
    sin::AbstractArray{<:AbstractFloat, 3},
    ls::AbstractArray{<:AbstractFloat, 3},
    hortot::Horizontal,
    state::State,
)
```

Apply the "horizontal part" of the linear operator (excluding the lower, center and upper diagonal) to the solution array `sin`.

Before the operator is applied, the boundary/halo values of `sin` are set, in the same way as in the method applying the total operator.

# Arguments

  - `sin`: Solution array.

  - `ls`: Result of applying the operator to the solution array.

  - `hortot`: Linear-operator mode.

  - `state`: Model state.

# See also

  - [`PinCFlow.Boundaries.set_zonal_boundaries_of_field!`](@ref)

  - [`PinCFlow.Boundaries.set_meridional_boundaries_of_field!`](@ref)

  - [`PinCFlow.MPIOperations.set_vertical_halos_of_field!`](@ref)
"""
function apply_operator! end

function apply_operator!(
    sin::AbstractArray{<:AbstractFloat, 3},
    ls::AbstractArray{<:AbstractFloat, 3},
    hortot::Total,
    state::State,
)
    (; namelists, domain) = state
    (; npz) = state.namelists.domain
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (;
        ac_b,
        al_b,
        ar_b,
        ab_b,
        af_b,
        ad_b,
        au_b,
        aru_b,
        ard_b,
        alu_b,
        ald_b,
        afu_b,
        afd_b,
        abu_b,
        abd_b,
        auu_b,
        add_b,
        aruu_b,
        ardd_b,
        aluu_b,
        aldd_b,
        afuu_b,
        afdd_b,
        abuu_b,
        abdd_b,
    ) = state.poisson.tensor
    (; s) = state.poisson.operator

    # Initialize auxiliary field.
    @ivy s[i0:i1, j0:j1, k0:k1] .= sin

    # Set boundaries of auxiliary field.
    if npz > 1
        set_zonal_boundaries_of_field!(s, namelists, domain; layers = (1, 1, 2))
        set_meridional_boundaries_of_field!(
            s,
            namelists,
            domain;
            layers = (1, 1, 2),
        )
        set_vertical_halos_of_field!(s, namelists, domain; layers = (1, 1, 2))
    else
        set_zonal_boundaries_of_field!(s, namelists, domain; layers = (1, 1, 1))
        set_meridional_boundaries_of_field!(
            s,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    #---------------------------------
    #         Loop over field
    #---------------------------------

    @ivy for k in 1:nz, j in 1:ny, i in 1:nx

        # Determine indices for s.
        is = i + i0 - 1
        js = j + j0 - 1
        ks = k + k0 - 1

        # ------------------ A(i+1,j,k) ------------------

        ar = ar_b[i, j, k]
        sr = s[is + 1, js, ks]

        # ------------------- A(i-1,j,k) --------------------

        al = al_b[i, j, k]
        sl = s[is - 1, js, ks]

        # -------------------- A(i,j+1,k) ----------------------

        af = af_b[i, j, k]
        sf = s[is, js + 1, ks]

        # --------------------- A(i,j-1,k) -----------------------

        ab = ab_b[i, j, k]
        sb = s[is, js - 1, ks]

        # --------------------- A(i,j,k+1) ------------------------

        au = au_b[i, j, k]
        su = s[is, js, ks + 1]

        # --------------------- A(i,j,k-1) ------------------------

        ad = ad_b[i, j, k]
        sd = s[is, js, ks - 1]

        # -------------------- A(i,j,k) --------------------------

        ac = ac_b[i, j, k]
        sc = s[is, js, ks]

        # -------------------- Apply operator ---------------------

        ls[i, j, k] =
            al * sl + ar * sr + af * sf + ab * sb + au * su + ad * sd + ac * sc

        # ----------------- A(i+1,j,k+1) -----------------

        aru = aru_b[i, j, k]
        sru = s[is + 1, js, ks + 1]

        # ----------------- A(i+1,j,k-1) -----------------

        ard = ard_b[i, j, k]
        srd = s[is + 1, js, ks - 1]

        # ----------------- A(i-1,j,k+1) -----------------

        alu = alu_b[i, j, k]
        slu = s[is - 1, js, ks + 1]

        # ----------------- A(i-1,j,k-1) -----------------

        ald = ald_b[i, j, k]
        sld = s[is - 1, js, ks - 1]

        # ----------------- A(i,j+1,k+1) -----------------

        afu = afu_b[i, j, k]
        sfu = s[is, js + 1, ks + 1]

        # ----------------- A(i,j+1,k-1) -----------------

        afd = afd_b[i, j, k]
        sfd = s[is, js + 1, ks - 1]

        # ----------------- A(i,j-1,k+1) -----------------

        abu = abu_b[i, j, k]
        sbu = s[is, js - 1, ks + 1]

        # ----------------- A(i,j-1,k-1) -----------------

        abd = abd_b[i, j, k]
        sbd = s[is, js - 1, ks - 1]

        # ------------------ A(i,j,k+2) -----------------

        auu = auu_b[i, j, k]
        suu = s[is, js, ks + 2]

        # ------------------ A(i,j,k-2) -----------------

        add = add_b[i, j, k]
        sdd = s[is, js, ks - 2]

        # ----------------- A(i+1,j,k+2) -----------------

        aruu = aruu_b[i, j, k]
        sruu = s[is + 1, js, ks + 2]

        # ----------------- A(i+1,j,k-2) -----------------

        ardd = ardd_b[i, j, k]
        srdd = s[is + 1, js, ks - 2]

        # ----------------- A(i-1,j,k+2) -----------------

        aluu = aluu_b[i, j, k]
        sluu = s[is - 1, js, ks + 2]

        # ----------------- A(i-1,j,k-2) -----------------

        aldd = aldd_b[i, j, k]
        sldd = s[is - 1, js, ks - 2]

        # ----------------- A(i,j+1,k+2) -----------------

        afuu = afuu_b[i, j, k]
        sfuu = s[is, js + 1, ks + 2]

        # ----------------- A(i,j+1,k-2) -----------------

        afdd = afdd_b[i, j, k]
        sfdd = s[is, js + 1, ks - 2]

        # ----------------- A(i,j-1,k+2) -----------------

        abuu = abuu_b[i, j, k]
        sbuu = s[is, js - 1, ks + 2]

        # ----------------- A(i,j-1,k-2) -----------------

        abdd = abdd_b[i, j, k]
        sbdd = s[is, js - 1, ks - 2]

        # Update operator.
        ls[i, j, k] +=
            aru * sru +
            ard * srd +
            alu * slu +
            ald * sld +
            afu * sfu +
            afd * sfd +
            abu * sbu +
            abd * sbd +
            auu * suu +
            add * sdd +
            aruu * sruu +
            ardd * srdd +
            aluu * sluu +
            aldd * sldd +
            afuu * sfuu +
            afdd * sfdd +
            abuu * sbuu +
            abdd * sbdd
    end
    return
end

function apply_operator!(
    sin::AbstractArray{<:AbstractFloat, 3},
    ls::AbstractArray{<:AbstractFloat, 3},
    hortot::Horizontal,
    state::State,
)
    (; namelists, domain) = state
    (; npz) = state.namelists.domain
    (; nx, ny, nz, i0, i1, j0, j1, k0, k1) = state.domain
    (;
        al_b,
        ar_b,
        ab_b,
        af_b,
        aru_b,
        ard_b,
        alu_b,
        ald_b,
        afu_b,
        afd_b,
        abu_b,
        abd_b,
        auu_b,
        add_b,
        aruu_b,
        ardd_b,
        aluu_b,
        aldd_b,
        afuu_b,
        afdd_b,
        abuu_b,
        abdd_b,
    ) = state.poisson.tensor
    (; s) = state.poisson.operator

    # Initialize auxiliary field.
    @ivy s[i0:i1, j0:j1, k0:k1] .= sin

    # Set boundaries of auxiliary field.
    if npz > 1
        set_zonal_boundaries_of_field!(s, namelists, domain; layers = (1, 1, 2))
        set_meridional_boundaries_of_field!(
            s,
            namelists,
            domain;
            layers = (1, 1, 2),
        )
        set_vertical_halos_of_field!(s, namelists, domain; layers = (1, 1, 2))
    else
        set_zonal_boundaries_of_field!(s, namelists, domain; layers = (1, 1, 1))
        set_meridional_boundaries_of_field!(
            s,
            namelists,
            domain;
            layers = (1, 1, 1),
        )
    end

    #---------------------------------
    #         Loop over field
    #---------------------------------

    @ivy for k in 1:nz, j in 1:ny, i in 1:nx

        # Determine indices for s.
        is = i + i0 - 1
        js = j + j0 - 1
        ks = k + k0 - 1

        # ------------------ A(i+1,j,k) ------------------

        ar = ar_b[i, j, k]
        sr = s[is + 1, js, ks]

        # ------------------- A(i-1,j,k) --------------------

        al = al_b[i, j, k]
        sl = s[is - 1, js, ks]

        # -------------------- A(i,j+1,k) ----------------------

        af = af_b[i, j, k]
        sf = s[is, js + 1, ks]

        # --------------------- A(i,j-1,k) -----------------------

        ab = ab_b[i, j, k]
        sb = s[is, js - 1, ks]

        ls[i, j, k] = al * sl + ar * sr + af * sf + ab * sb

        # ----------------- A(i+1,j,k+1) -----------------

        aru = aru_b[i, j, k]
        sru = s[is + 1, js, ks + 1]

        # ----------------- A(i+1,j,k-1) -----------------

        ard = ard_b[i, j, k]
        srd = s[is + 1, js, ks - 1]

        # ----------------- A(i-1,j,k+1) -----------------

        alu = alu_b[i, j, k]
        slu = s[is - 1, js, ks + 1]

        # ----------------- A(i-1,j,k-1) -----------------

        ald = ald_b[i, j, k]
        sld = s[is - 1, js, ks - 1]

        # ----------------- A(i,j+1,k+1) -----------------

        afu = afu_b[i, j, k]
        sfu = s[is, js + 1, ks + 1]

        # ----------------- A(i,j+1,k-1) -----------------

        afd = afd_b[i, j, k]
        sfd = s[is, js + 1, ks - 1]

        # ----------------- A(i,j-1,k+1) -----------------

        abu = abu_b[i, j, k]
        sbu = s[is, js - 1, ks + 1]

        # ----------------- A(i,j-1,k-1) -----------------

        abd = abd_b[i, j, k]
        sbd = s[is, js - 1, ks - 1]

        # ------------------ A(i,j,k+2) -----------------

        auu = auu_b[i, j, k]
        suu = s[is, js, ks + 2]

        # ------------------ A(i,j,k-2) -----------------

        add = add_b[i, j, k]
        sdd = s[is, js, ks - 2]

        # ----------------- A(i+1,j,k+2) -----------------

        aruu = aruu_b[i, j, k]
        sruu = s[is + 1, js, ks + 2]

        # ----------------- A(i+1,j,k-2) -----------------

        ardd = ardd_b[i, j, k]
        srdd = s[is + 1, js, ks - 2]

        # ----------------- A(i-1,j,k+2) -----------------

        aluu = aluu_b[i, j, k]
        sluu = s[is - 1, js, ks + 2]

        # ----------------- A(i-1,j,k-2) -----------------

        aldd = aldd_b[i, j, k]
        sldd = s[is - 1, js, ks - 2]

        # ----------------- A(i,j+1,k+2) -----------------

        afuu = afuu_b[i, j, k]
        sfuu = s[is, js + 1, ks + 2]

        # ----------------- A(i,j+1,k-2) -----------------

        afdd = afdd_b[i, j, k]
        sfdd = s[is, js + 1, ks - 2]

        # ----------------- A(i,j-1,k+2) -----------------

        abuu = abuu_b[i, j, k]
        sbuu = s[is, js - 1, ks + 2]

        # ----------------- A(i,j-1,k-2) -----------------

        abdd = abdd_b[i, j, k]
        sbdd = s[is, js - 1, ks - 2]

        # Update operator.
        ls[i, j, k] +=
            aru * sru +
            ard * srd +
            alu * slu +
            ald * sld +
            afu * sfu +
            afd * sfd +
            abu * sbu +
            abd * sbd +
            auu * suu +
            add * sdd +
            aruu * sruu +
            ardd * srdd +
            aluu * sluu +
            aldd * sldd +
            afuu * sfuu +
            afdd * sfdd +
            abuu * sbuu +
            abdd * sbdd
    end
    return
end
