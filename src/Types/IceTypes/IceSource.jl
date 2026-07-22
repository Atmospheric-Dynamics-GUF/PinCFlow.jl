struct IceSource{A <: AbstractArray{<:AbstractFloat, 3}}
    n_homsource :: A
    q_homsource :: A
    qvsource :: A
    n_insource :: A
    n_hetsource :: A
    q_hetsource :: A
end

function IceSource(namelists::Namelists, domain::Domain)    
    (; ice_setup, nucleation) = namelists.ice

    return IceSource(domain, ice_setup, nucleation)
end

# no ice
function IceSource(domain::Domain, ice_setup::NoIce, nucleation::AbstractNucleation)
    
    n_homsource = zeros(0, 0, 0)
    q_homsource = zeros(0, 0, 0)
    qvsource = zeros(0, 0, 0)
    n_insource = zeros(0, 0, 0)
    n_hetsource = zeros(0, 0, 0)
    q_hetsource = zeros(0, 0, 0)

    return IceSource(n_homsource, q_homsource, qvsource, n_insource, n_hetsource, q_hetsource)
end

# both nucleation
function IceSource(domain::Domain, ice_setup::AbstractIce, nucleation::BothNucleation)
    (; nxx, nyy, nzz) = domain

    n_homsource = zeros(nxx, nyy, nzz)
    q_homsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)
    n_insource = zeros(nxx, nyy, nzz)
    n_hetsource = zeros(nxx, nyy, nzz)
    q_hetsource = zeros(nxx, nyy, nzz)

    return IceSource(n_homsource, q_homsource, qvsource, n_insource, n_hetsource, q_hetsource)
end

# only homogeneous nucleation
function IceSource(domain::Domain, ice_setup::AbstractIce, nucleation::HomogeneousOnly)
    (; nxx, nyy, nzz) = domain
    #=
    n_homsource = zeros(nxx, nyy, nzz)
    q_homsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)
    n_insource = zeros(nxx, nyy, nzz)
    n_hetsource = zeros(nxx, nyy, nzz)
    q_hetsource = zeros(nxx, nyy, nzz)
    =#
    n_homsource = zeros(nxx, nyy, nzz)
    q_homsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)
    n_insource = zeros(0, 0, 0)
    n_hetsource = zeros(0, 0, 0)
    q_hetsource = zeros(0, 0, 0)
    

    return IceSource(n_homsource, q_homsource, qvsource, n_insource, n_hetsource, q_hetsource)
end

# only heterogeneous nucleation
function IceSource(domain::Domain, ice_setup::AbstractIce, nucleation::HeterogeneousOnly)
    (; nxx, nyy, nzz) = domain
    #=
    n_homsource = zeros(nxx, nyy, nzz)
    q_homsource = zeros(nxx, nyy, nzz)
    qvsource = zeros(nxx, nyy, nzz)
    n_insource = zeros(nxx, nyy, nzz)
    n_hetsource = zeros(nxx, nyy, nzz)
    q_hetsource = zeros(nxx, nyy, nzz)
    =#
    n_homsource = zeros(0, 0, 0)
    q_homsource = zeros(0, 0, 0)
    qvsource = zeros(nxx, nyy, nzz)
    n_insource = zeros(nxx, nyy, nzz)
    n_hetsource = zeros(nxx, nyy, nzz)
    q_hetsource = zeros(nxx, nyy, nzz)
    

    return IceSource(n_homsource, q_homsource, qvsource, n_insource, n_hetsource, q_hetsource)
end
