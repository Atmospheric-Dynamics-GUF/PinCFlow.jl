function compute_msgwam_ice! end

function compute_msgwam_ice!(state::State)
	(; wkb_mode) = state.namelists.wkb
	compute_msgwam_ice!(state, wkb_mode)
	return
end

function compute_msgwam_ice!(state::State, wkb_mode::NoWKB)
	return
end

function compute_msgwam_ice!(state::State, wkb_mode::Union{SteadyState, SingleColumn, MultiColumn})
	ice_setup = state.namelists.ice.ice_setup	
	compute_msgwam_ice!(state, ice_setup)
	return
end

function compute_msgwam_ice!(state::State, ice_setup::NoIce)
	return
end