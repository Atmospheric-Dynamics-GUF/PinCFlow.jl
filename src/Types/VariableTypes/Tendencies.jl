struct Tendencies{
    A <: AbstractArray{<:AbstractFloat, 3},
    B <: AbstractArray{<:AbstractFloat, 3},
}
    drho::A
    drhop::A
    du::A
    dv::A
    dw::A
    dpip::A
    dp::B
end

function Tendencies(namelists::Namelists, domain::Domain)
    (; model) = namelists.setting
    return Tendencies(domain, model)
end

function Tendencies(domain::Domain, model::AbstractModel)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip) = (zeros(nxx, nyy, nzz) for i in 1:6)
    dp = zeros(0, 0, 0)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end

function Tendencies(domain::Domain, model::Compressible)
    (; nxx, nyy, nzz) = domain

    # Initialize the tendencies.
    (drho, drhop, du, dv, dw, dpip, dp) = (zeros(nxx, nyy, nzz) for i in 1:7)

    # Return a Variables instance.
    return Tendencies(drho, drhop, du, dv, dw, dpip, dp)
end
