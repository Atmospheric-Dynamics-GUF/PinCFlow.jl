struct Backups{A <: AbstractArray{<:AbstractFloat, 3}}
    rhoold::A
    rhopold::A
    uold::A
    vold::A
    wold::A
end

function Backups(domain::Domain)

    # Get parameters.
    (; nxx, nyy, nzz) = domain

    # Initialize the backups.
    (rhoold, rhopold, uold, vold, wold) = (zeros(nxx, nyy, nzz) for i in 1:5)

    # Return a Backups instance.
    return Backups(rhoold, rhopold, uold, vold, wold)
end
