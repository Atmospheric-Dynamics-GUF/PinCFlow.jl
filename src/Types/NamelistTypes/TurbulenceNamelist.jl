struct TurbulenceNamelist{A <: AbstractTurbulence}
    turbulencesetup::A
end

function TurbulenceNamelist(; turbulencesetup = NoTurbulence())
    return TurbulenceNamelist(turbulencesetup)
end
