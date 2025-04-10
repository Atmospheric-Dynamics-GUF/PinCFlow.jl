struct Tendencies{A <: AbstractArray{<:AbstractFloat, 3}}
    drho::A
    drhop::A
    du::A
    dv::A
    dw::A
    dpip::A
end

function Tendencies(domain::Domain)

    # Get parameters.
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip) = (zeros(nxx, nyy, nzz) for i in 1:6)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip)
end
