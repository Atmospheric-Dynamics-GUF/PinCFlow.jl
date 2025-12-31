function compute_m1m2 end


function compute_m1m2(
    kpr::AbstractFloat, 
    kpr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr::AbstractFloat,
    res_type::Sum,
    sol_branch::Sum
    )::NTuple{2, <:AbstractFloat}

    f1 = f2 = m1 = m2 = 0.0

    f1 =  mr / (2 * kpr)
    f2 = kpr + kpr1 + kpr2
    m1 = f1 * (f2 + sqrt(f2^2 - 4 * kpr * kpr1))
    m2 = mr - m1

    return (m1, m2)

end


function compute_m1m2(
    kpr::AbstractFloat, 
    kpr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr::AbstractFloat,
    res_type::Sum,
    sol_branch::Difference
    )::NTuple{2, <:AbstractFloat}
    
    f1 = f2 = m1 = m2 = 0.0

    f1 =  mr / (2 * kpr)
    f2 = kpr - kpr1 - kpr2
    m1 = f1 * (f2 - sqrt(f2^2 + 4 * kpr * kpr1))
    m2 = mr - m1

    return (m1, m2)

end

function compute_m1m2(
    kpr::AbstractFloat, 
    kpr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr::AbstractFloat,
    res_type::Difference,
    sol_branch::Sum
    )::NTuple{2, <:AbstractFloat}

    f1 = f2 = m1 = m2 = 0.0
    
    f1 =  mr / (2 * kpr)
    f2 = kpr + kpr1 + kpr2
    m1 = f1 * (f2 - sqrt(f2^2 - 4 * kpr * kpr1))
    m2 = m1 - mr

    return (m1, m2)

end

function compute_m1m2(
    kpr::AbstractFloat, 
    kpr1::AbstractFloat,
    kpr2::AbstractFloat,
    mr::AbstractFloat,
    res_type::Difference,
    sol_branch::Difference
    )::NTuple{2, <:AbstractFloat}

    f1 = f2 = m1 = m2 = 0.0
    
    f1 =  mr / (2 * kpr)
    f2 = kpr - kpr1 + kpr2
    m1 = f1 * (f2 - sqrt((-f2)^2 + 4 * kpr * kpr1))
    m2 = m1 - mr

    return (m1, m2)

end