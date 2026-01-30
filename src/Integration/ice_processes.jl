# contains additional ice process functions 
#
# processes are not physically accurate
# and are only intended to produce data
# for training and testing ML models

# z_factor to localize source/sink terms in vertical
function z_factor(constants::Constants, grid::Grid, i::Int, j::Int, k::Int)
    (; zc) = grid
    (; lref) = constants
    
    # **********************
	# ISSRegion as IC
	# **********************    

	# center ISSR
	z0_issr = 8.e3 # [m]
	# vertical width ISSR (standard deviation of gaussian dist.)
	sig_issr = 4.e3 # [m]
	# max value S in ISSR
	S_issr = 1.45 #iceconstants.S_c 

	#nondim.
	z0_issr = z0_issr / lref
	sig_issr = sig_issr / lref

	#define upper/lower bounds of ISSR
	zMin_issr = z0_issr - sig_issr
	zMax_issr = z0_issr + sig_issr

    return 0.5 * (tanh( (zc[i, j, k] - zMin_issr) / (0.1 * sig_issr) ) - tanh( (zc[i, j, k] - zMax_issr) / (0.1 * sig_issr) ))
end

# sedimentation sink term for q
function q_sedimentation(q::AbstractFloat, n::AbstractFloat, namelists::Namelists)
    (; tau) = namelists.ice.tau_q_sink
    return - 1.0 / tau * q^(5.0 / 3.0) * n^(-2.0 / 3.0)
end

# vapor source term for qv
function qv_forcing(qv::AbstractFloat, q::AbstractFloat, namelists::Namelists, constants::Constants, grid::Grid, i::Int, j::Int, k::Int)
    (; tau_qv_source) = namelists.ice.tau_qv_source

    # equilibrium vapor mixing ratio
    qref = 1.0 # reference mixing ratio
    qv_eq = 0.052 # equilibrium vapor mixing ratio

    return (qv_eq - qv) / tau_qv_source * exp(-q / qref) * z_factor(constants, grid, i, j, k)
end