struct ResManifold{A <: AbstractVector{<:AbstractFloat}}
    k1s::A  
    k2s::A
    k1d::A  
    k2d::A
    m1sp::A  
    m2sp::A
    m1dp::A 
    m2dp::A
    m1sn::A
    m2sn::A
    m1dn::A
    m2dn::A
end

function ResManifold(wkb_mode::Union{NoWKB, SteadyState, SingleColumn, MultiColumn},
    triad_mode::NoTriad)::ResManifold

    return ResManifold(zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0), zeros(0))

end


function ResManifold(l_sum::Integer,
   l_diff::Integer,
   wkb_mode::Union{SteadyState, SingleColumn, MultiColumn},
   triad_mode::Union{Triad2D, Triad3DIso})::ResManifold

   ni_sum = Int(2 * l_sum - 1)
    k1s = zeros(ni_sum) 
    k2s = zeros(ni_sum)
    k1d = zeros(l_diff)  
    k2d = zeros(l_diff)
    m1sp = zeros(ni_sum)
    m2sp = zeros(ni_sum)
    m1dp = zeros(l_diff) 
    m2dp = zeros(l_diff)
    m1sn = zeros(ni_sum)
    m2sn = zeros(ni_sum)
    m1dn = zeros(l_diff) 
    m2dn = zeros(l_diff)

   
   return ResManifold(k1s, k2s, k1d, k2d, m1sp, m2sp, m1dp, m2dp, m1sn, m2sn, m1dn, m2dn)
    
end