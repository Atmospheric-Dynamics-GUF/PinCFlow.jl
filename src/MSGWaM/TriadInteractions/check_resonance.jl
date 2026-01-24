function check_resonance end

function check_resonance(
    kpr::AbstractFloat, 
    mr::AbstractFloat,
    kpr1::AbstractFloat,
    mr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr2::AbstractFloat,
    res_type::Sum,
    )

    reso_offening = abs(abs(kpr1) / abs(mr1) + abs(kpr2) / abs(mr2) - abs(kpr) / abs(mr))
    if reso_offening > 1.0E-5
        error("Error in resonance, resonance offening for the sum interaction is", reso_offening)
    end
    return 

end

function check_resonance(
    kpr::AbstractFloat, 
    mr::AbstractFloat,
    kpr1::AbstractFloat,
    mr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr2::AbstractFloat,
    res_type::Difference,
    )

    reso_offening = abs(kpr1) / abs(mr1) - abs(kpr2) / abs(mr2) - abs(kpr) / abs(mr) 
    if reso_offening > 1.0E-5
        error("Error in resonance, resonance offening for the difference interaction is", reso_offening)
    end
    return 

end