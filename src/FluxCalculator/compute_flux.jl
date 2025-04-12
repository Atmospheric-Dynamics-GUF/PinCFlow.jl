function compute_flux(
    usurf::AbstractFloat,
    phiup::AbstractFloat,
    phidown::AbstractFloat,
)
    if usurf > 0.0
        return usurf * phiup
    else
        return usurf * phidown
    end
end
