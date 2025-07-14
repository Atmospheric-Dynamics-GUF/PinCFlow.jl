struct IceAuxiliaries{A <: AbstractArray{<:AbstractFloat, 3}}
    initialn::A
    initialq::A
    initialqv::A
end

function IceAuxiliaries(icepredictands::IcePredictands)
    initialn = copy(getfield(icepredictands, :n))
    initialq = copy(getfield(icepredictands, :q))
    initialqv = copy(getfield(icepredictands, :qv))

    return IceAuxiliaries(initialn, initialq, initialqv)
end

struct IcePhysics{A <: AbstractFloat}       
     J_nuc :: A
     B_nuc :: A
end 
function IcePhysics()
    return(J_nuc=999., B_nuc=999.)
end

# """

# struct IcePhysics{A <: AbstractFloat}
#       J_nuc, B_nuc, epsil0, meanMassIce, mref :: A
     
#       PsatIceRef, Dep, Dep0, DepS :: A
#       L_ice, R_v, thetaRef_trp :: A
#       dt_ice :: A

#       thetaRefRatio :: A 
#       L_hat, Li_hat, epsil0hat :: A

# end    

# function IcePhysics(J) 
      
#     return IcePhysics(
#         J_nuc = 999., 
#         B_nuc = 999.,
#         epsil0 = 0.62, #\approx Mole_mass_water/Mole_mass_dryAir
#         meanMassIce = 1.E-12, #mean mass ice crystals [kg]
#         mRef = 999., # reference mass
#         PsatIceRef = 1., #reference saturation pressure [Pa]
#         Dep = 999.,
#         Dep0 = 999.,
#         DepS = 999., # deposition coefficient
#         L_ice = 2.8E6, # constant latent heat  ice [J/kg]
#         R_v  = 461., # specific gas constant for water vapor [J/kg/K]
#         thetaRef_trp = 210., # reference temperature in the tropopause regio[K]
#         dt_ice = 0.1, # length of the microphysical time step [s]
#         thetaRefRatio = 999.,
#         L_hat =  999.,
#         Li_hat   = 999.,
#         epsil0hat     = 999.,
#     )
# end   
# """
# """
# function initialize_ice(
#     state::State,
#     icesetup::AbstractIce,
# )
#     (; icepredictands, iceauxiliaries) = state.ice

#     # Initialize ice predictands.
#     #initialize_ice_predictands!(icepredictands, icesetup)

#     # Initialize ice auxiliaries.
#     #initialize_ice_auxiliaries!(iceauxiliaries, icepredictands)

#     return
# end

# function explicit_integration_rhs_ice!(state::State, icesetup::NoIce)
#     # No ice setup, no ice integration.
#     return

# end    

# function explicit_integration_rhs_ice!(
#     state::State,
#     icesetup::AbstractIce,
# )
#     (; icepredictands, icetendencies, iceauxiliaries) = state.ice

#     # Compute ice tendencies.
#     #compute_ice_tendencies!(icetendencies, icepredictands, icesetup)

#     # Update ice predictands.
#     #update_ice_predictands!(icepredictands, icetendencies)

#     # Update ice auxiliaries.
#     #update_ice_auxiliaries!(iceauxiliaries, icepredictands)

#     return
# end
# """