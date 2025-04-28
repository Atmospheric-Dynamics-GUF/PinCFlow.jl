struct Auxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    phi::A
end

function Auxiliaries(domain::Domain)
    (; nxx, nyy, nzz) = domain

    # Initialize the auxiliaries.
    phi = zeros(nxx, nyy, nzz)

    # Return an Auxiliaries instance.
    return Auxiliaries(phi)
end
