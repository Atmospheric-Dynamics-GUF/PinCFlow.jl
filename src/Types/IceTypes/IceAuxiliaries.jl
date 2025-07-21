struct IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    iaux1::A
    iaux2::A
    iaux3::A
end

function IceAuxiliaries(icepredictands::IcePredictands)
    iaux1 = copy(getfield(icepredictands, :n))
    iaux2 = copy(getfield(icepredictands, :q))
    iaux3 = copy(getfield(icepredictands, :qv))

    return IceAuxiliaries(iaux1, iaux2, iaux3)
end
