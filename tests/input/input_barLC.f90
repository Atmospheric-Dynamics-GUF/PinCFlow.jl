!------------------------------------------------------------------------------!
!                              PincFlow Namelists                              !
!------------------------------------------------------------------------------!

&domain
sizeX                    =                      21 ! Cells in x
sizeY                    =                      84 ! Cells in y
sizeZ                    =                      75 ! Cells in z
nbx                      =                       2 ! Halo/ghost cells in x
nby                      =                       2 ! Halo/ghost cells in y
nbz                      =                       2 ! Ghost cells in z
lx_dim(0)                =                  0.0E+0 ! Minimum of x
lx_dim(1)                =                  4.2E+6 ! Maximum of x
ly_dim(0)                =                 -8.4E+6 ! Minimum of y
ly_dim(1)                =                  8.4E+6 ! Maximum of y
lz_dim(0)                =                  0.0E+0 ! Minimum of z
lz_dim(1)                =                  1.5E+5 ! Maximum of z
nprocx                   =                {nprocx} ! Processors in x
nprocy                   =                {nprocy} ! Processors in y
&end

&variables
include_ice              =                       F ! Ice microphysics
                                                   ! parameterization
include_tracer           =                       F ! Tracer equation
include_testoutput       =                       F ! ...
&end

&outputList
runName                  =                 'barLC' ! run name for netCDF file
atmvarOut(1)             =                   'rho' ! Atmospheric output
                                                   ! variable
atmvarOut(2)             =                     'u' ! Atmospheric output
                                                   ! variable
atmvarOut(3)             =                     'v' ! Atmospheric output
                                                   ! variable
atmvarOut(4)             =                     'w' ! Atmospheric output
                                                   ! variable
atmvarOut(5)             =                    'pi' ! Atmospheric output
                                                   ! variable
atmvarOut(6)             =                 'theta' ! Atmospheric output
                                                   ! variable
outputType               =                  'time' ! 'timeStep' or 'time'
nOutput                  =                       1 ! Output every nOutput
                                                   ! time steps for
                                                   ! outputType = 'timeStep'
maxIter                  =                     100 ! Stop after maxIter time
                                                   ! steps for outputType =
                                                   ! 'timeStep'
outputTimeDiff           =                1.728E+5 ! Output every
                                                   ! outputTimeDiff seconds
                                                   ! for outputType = 'time'
maxTime                  =               1.2096E+6 ! Stop after maxTime
                                                   ! seconds for outputType =
                                                   ! 'time'
iIn                      =                      97 ! restart at time-step iIn
restart                  =                       F ! Restart the model from
                                                   ! state in previous
                                                   ! simulation
restartIN(1)             =                   'rho' ! Restart input variable
restartIN(2)             =                    'us' ! Restart input variable
restartIN(3)             =                    'vs' ! Restart input variable
restartIN(4)             =                    'ws' ! Restart input variable
detailedinfo             =                       F ! Provide info on the
                                                   ! final state of Poisson
                                                   ! solver
RHS_diagnostics          =                       T ! Provide info about the
                                                   ! RHS of Poisson equation
fancy_namelists          =                       T ! Write all namelists with
                                                   ! comments
&end

&debuggingList
verbose                  =                       F ! ...
dtMin_dim                =                  1.0E-5 ! Stop if dt < dtMin
&end

&testCaseList
testCase                 =         'baroclinic_LC' ! Predefined test case
&end

&monochromeWave
lambda_dim               =                  6.0E+3 ! Vertical wavelength
&end

&wavePacket
wavePacketType           =                       1 ! 1 = Gaussian, 2 = Cosine
wavePacketDim            =                       2 ! 1 = 1D, 2 = 2D, 3 = 3D
lambdaX_dim              =                  3.0E+5 ! Wavelength in x (0.0 for
                                                   ! infinite wavelength)
lambdaY_dim              =                  0.0E+0 ! Wavelength in y (0.0 for
                                                   ! infinite wavelength)
lambdaZ_dim              =                 -1.0E+3 ! Wavelength in z
amplitudeFactor          =                  0.0E+0 ! Normalized buoyancy
                                                   ! amplitude
x0_dim                   =                  4.5E+5 ! Wave-packet center in x
y0_dim                   =                  2.0E+4 ! Wave-packet center in y
z0_dim                   =                  3.0E+4 ! Wave-packet center in z
sigma_dim                =                  5.0E+3 ! Vertical width of
                                                   ! Gaussian wave packet
sigma_hor_dim            =                  1.5E+6 ! Cosine distribution
                                                   ! width in x (0.0 for
                                                   ! infinite width)
amp_mod_x                =                  1.0E+0 ! Fractional amplitude
                                                   ! modulation in x
sigma_hor_yyy_dim        =                  0.0E+0 ! Cosine distribution
                                                   ! width in y (0.0 for
                                                   ! infinite width)
amp_mod_y                =                  1.0E+0 ! Fractional amplitude
                                                   ! modulation in y
L_cos_dim                =                  1.0E+4 ! Half width of vertical
                                                   ! cosine profile of wave
                                                   ! packet
omiSign                  =                       1 ! Frequency branch
u0_jet_dim               =                  0.0E+0 ! Maximum amplitude of jet
z0_jet_dim               =                  5.0E+4 ! Center of jet
L_jet_dim                =                  5.0E+3 ! Half width of vertical
                                                   ! cosine profile of jet
inducedwind              =                       F ! ...
&end

&LagrangeRayTracing
xrmin_dim                =                  0.0E+0 ! Left bound of initial
                                                   ! rays
xrmax_dim                =                  5.0E+6 ! Right bound of initial
                                                   ! rays
yrmin_dim                =                  0.0E+0 ! Backward bound of
                                                   ! initial rays
yrmax_dim                =                  5.0E+4 ! Forward bound of initial
                                                   ! rays
zrmin_dim                =                  3.0E+3 ! Bottom bound of initial
                                                   ! rays
zrmax_dim                =                  7.0E+4 ! Top bound of initial rays
nrxl                     =                       1 ! Initial ray volumes
                                                   ! within dx
nryl                     =                       1 ! Initial ray volumes
                                                   ! within dy
nrzl                     =                       1 ! Initial ray volumes
                                                   ! within dz
fac_dk_init              =                  1.0E-1 ! Initial fraction dk/k
fac_dl_init              =                  1.0E-1 ! Initial fraction dl/l
fac_dm_init              =                  1.0E-1 ! Initial fraction dm/m
nrk_init                 =                       2 ! Initial ray volumes
                                                   ! within dk
nrl_init                 =                       1 ! Initial ray volumes
                                                   ! within dl
nrm_init                 =                       2 ! Initial ray volumes
                                                   ! within dm
nsmth_wkb                =                       2 ! Width of smoothing
                                                   ! operator for mean flow
                                                   ! tendencies
lsmth_wkb                =                       T ! Smoothing operator for
                                                   ! mean flow tendencies
sm_filter                =                       2 ! 1 = Box filter, 2 =
                                                   ! Shapiro filter
lsaturation              =                       T ! Switch for saturation
                                                   ! scheme
alpha_sat                =                  1.0E+0 ! Saturation threshold
steady_state             =                       F ! Steady-state mode
case_wkb                 =                       3 ! 1 = Gaussian wave
                                                   ! packet, 2 = cosine wave
                                                   ! packet, 3 = mountain
                                                   ! wave
amp_wkb                  =                  5.0E-1 ! Relative amplitude of
                                                   ! wave packet
wlrx_init                =                  1.0E+5 ! Wavelength in x of wave
                                                   ! packet
wlry_init                =                  0.0E+0 ! Wavelength in y of wave
                                                   ! packet
wlrz_init                =                  1.0E+3 ! Wavelength in z of wave
                                                   ! packet
xr0_dim                  =                  2.5E+6 ! Wave packet center in x
yr0_dim                  =                  2.5E+4 ! Wave packet center in y
zr0_dim                  =                  3.0E+4 ! Wave packet center in z
sigwpx_dim               =                  1.0E+6 ! Wave packet width in x
                                                   ! (0.0 for infinite width)
sigwpy_dim               =                  0.0E+0 ! Wave packet width in y
                                                   ! (0.0 for infinite width)
sigwpz_dim               =                  5.0E+3 ! Wave packet width in z
                                                   ! (0.0 for infinite width)
branchr                  =                       1 ! Frequency branch
lindUinit                =                       F ! Induced wind at initial
                                                   ! time
topographyTime_wkb       =                  0.0E+0 ! WKB topography growth
                                                   ! time
mountainHeight_wkb_dim   =                  5.0E+2 ! WKB mountain height
mountainWidth_wkb_dim    =                  1.0E+6 ! WKB mountain half width
mountain_case_wkb        =                       1 ! WKB topography shape
range_factor_wkb         =                  1.0E+1 ! Ratio between large and
                                                   ! small WKB scales
blocking                 =                       F ! Simple blocked-layer
                                                   ! scheme
nwm                      =                       1 ! Number of initial wave
                                                   ! modes
launch_level             =                       0 ! Ray-volume launch level
launch_algorithm         =                  'clip' ! Ray-volume launch
                                                   ! algorithm
zmin_wkb_dim             =                  0.0E+0 ! Minimum altitude for
                                                   ! wave-mean-flow
                                                   ! interaction
nray_fac                 =                       5 ! Maximum multiplication
                                                   ! factor per spectral
                                                   ! dimension
cons_merge               =                    'en' ! Conserved quantity in
                                                   ! ray-volume merging ('wa'
                                                   ! = wave action, 'en' =
                                                   ! wave energy)
nRayOutput               =                      10 ! Number of dominant ray
                                                   ! volumes in output
&end

&bubble
dTheta0_dim              =                  6.0E+0 ! ...
xRadius_dim              =                  2.5E+3 ! ...
zRadius_dim              =                  2.5E+3 ! ...
xCenter_dim              =                  0.0E+0 ! ...
yCenter_dim              =                  0.0E+0 ! ...
zCenter_dim              =                  0.0E+0 ! ...
zExcentricity            =                  1.0E+0 ! ...
rhoCenter_dim            =                  1.0E+3 ! ...
&end

&robert_bubble
dTheta1_dim              =                  5.0E-1 ! Potential temperature
                                                   ! offset
a1_dim                   =                  1.5E+2 ! Radius of plateau
sigma1_dim               =                  5.0E+1 ! Gaussian edge profile
xCenter1_dim             =                  0.0E+0 ! ...
zCenter1_dim             =                  3.0E+2 ! ...
dTheta2_dim              =                 -1.5E-1 ! ...
a2_dim                   =                  0.0E+0 ! ...
sigma2_dim               =                  5.0E+1 ! ...
xCenter2_dim             =                  6.0E+1 ! ...
zCenter2_dim             =                  6.4E+2 ! ...
&end

&mountainwavelist
u_relax                  =                  7.5E+1 ! Zonal relaxation wind
v_relax                  =                  0.0E+0 ! Meridional relaxation
                                                   ! wind
w_relax                  =                  0.0E+0 ! Vertical relaxation wind
t_relax                  =                 8.64E+4 ! Relaxation time
t_ramp                   =                  3.6E+3 ! Not used at the moment
xextent_relax            =                  3.0E+6 ! Zonal extent of
                                                   ! relaxation region
yextent_relax            =                  0.0E+0 ! Meridional extent of
                                                   ! relaxation region
wind_relaxation          =                       F ! Relaxation switch
surface_layer_depth      =                  0.0E+0 ! Surface-layer depth
&end

&baroclinic_LC
zero_initial_state       =                       F ! ...
z_trpp0_dim              =                  8.0E+3 ! Mean tropopause height
z_baro_dim               =                  2.4E+4 ! Altitude above which the
                                                   ! atmosphere is barotropic
thet0_dim                =                 2.73E+2 ! Characteristic potential
                                                   ! temperature
ntrp_dim                 =                  1.0E-2 ! Buoyancy frequency of
                                                   ! the troposphere
nstr_dim                 =                 2.45E-2 ! Buoyancy frequency of
                                                   ! the stratosphere
jwdth_dim                =                  2.1E+6 ! Jet half width
kaptpp                   =                  7.0E+1 ! Jet slope
add_ptptb                =                       T ! Add local potential-
                                                   ! temperature perturbation
ptptb_x_dim              =                  2.1E+6 ! Location in x of local
                                                   ! potential-temperature
                                                   ! perturbation
ptptb_y_dim              =                  4.2E+6 ! Location in y of local
                                                   ! potential-temperature
                                                   ! perturbation
ptptb_z_dim              =                  1.1E+4 ! Location in z of local
                                                   ! potential-temperature
                                                   ! perturbation
ptptb_dh_dim             =                  1.0E+6 ! Horizontal width of
                                                   ! local potential-
                                                   ! temperature perturbation
ptptb_dz_dim             =                  4.0E+3 ! Vertical width of local
                                                   ! potential-temperature
                                                   ! perturbation
ptptb_amp_dim            =                  3.0E-1 ! Amplitude of local
                                                   ! potential-temperature
                                                   ! perturbation
add_noise                =                       F ! Add noise to the initial
                                                   ! potential temperature
proc_noise               =                  1.0E-3 ! Relative amplitude of
                                                   ! potential-temperature
                                                   ! noise
tau_relax                =                 8.64E+5 ! ...
tau_relax_low            =                 8.64E+5 ! ...
sigma_tau                =                  6.0E-2 ! Relative thickness of
                                                   ! relaxation profile
tau_jet                  =                  7.2E+3 ! Time for jet formation
Sponge_Rel_Bal_Type      =                   'env' ! Relaxation state: 'hyd'
                                                   ! = hydrostatic balance,
                                                   ! 'env' = geostrophic
                                                   ! balance
ta_hs_dim                =                3.456E+6 ! Thermal-relaxation time
                                                   ! scale outside tropical
                                                   ! boundary layer (zero
                                                   ! means infinity)
ts_hs_dim                =                3.456E+5 ! Thermal-relaxation time
                                                   ! scale in tropical
                                                   ! boundary layer (zero
                                                   ! means infinity)
tf_hs_dim                =                 8.64E+4 ! Boundary-layer Rayleigh-
                                                   ! damping  time scale
sigb_hs                  =                  7.0E-1 ! Sigma of boundary-layer
                                                   ! top
&end

&baroclinic_ID
bar_sigma_y              =                  1.0E-1 ! ...
u_strength               =                  3.0E+1 ! Jet strength
dTh_atm                  =                  2.0E+1 ! Potential-temperature
                                                   ! difference between poles
                                                   ! and tropics
init_2Dto3D              =                       F ! Initialize 3D state from
                                                   ! 2D balanced state
init_bal                 =             'geostr_id' ! ...
lastrecordnum            =                      88 ! Last record in 2D output
fileinitstate2D          =    'pf_all_2D3D_in.dat' ! File name with 2D
                                                   ! balanced state
output_theta_bgr         =                       F ! Output environmental
                                                   ! potential temperature
output_br_vais_sq        =                       T ! Output environmental
                                                   ! squared buoyancy
                                                   ! frequency
output_heat              =                       F ! ...
balance_eq               =                 'PI_an' ! ...
output_rho_bgr           =                       T ! Output environmental
                                                   ! density
&end

&modelList
model                    = 'pseudo_incompressible' ! Dynamic equations
vert_theta               =                  9.0E+1 ! Rotation about x
vert_alpha               =                  0.0E+0 ! Rotation about y
&end

&solverList
cfl                      =                1.667E-1 ! CFL number
cfl_wave                 =                  5.0E-1 ! WKB CFL number
dtMax_dim                =                  1.8E+3 ! Maximum time step
tStepChoice              =                   'cfl' ! 'fix' or 'cfl'
timeScheme               =          'semiimplicit' ! 'LS_Will_RK3' or
                                                   ! 'semiimplicit'
auxil_equ                =                       F ! Buoyancy equation
fluxType                 =                'upwind' ! 'ILES', 'central' or
                                                   ! 'upwind'
reconstType              =                 'MUSCL' ! 'MUSCL', 'constant',
                                                   ! 'SALD' or 'ALDM'
musclType                =                'muscl1' ! 'muscl1' or 'muscl2'
limiterType1             =             'MCVariant' ! 'minmod', 'MCVariant' or
                                                   ! 'Cada'
TurbScheme               =                       F ! Turbulence scheme
turb_dts                 =                  5.0E+3 ! Turbulent damping time
DySmaScheme              =                       F ! Dynamic Smagorinsky
                                                   ! scheme
dtWave_on                =                       T ! Limit time step by
                                                   ! inverse buoyancy
                                                   ! frequency
heatingONK14             =                       T ! Heating as in O'Neill
                                                   ! and Klein (2014)
dens_relax               =                       F ! Heating by density
                                                   ! relaxation
shap_dts_fac             =                  1.0E+1 ! Shapiro-filter damping
                                                   ! time
n_shap                   =                       4 ! Order of Shapiro filter
&end

&poissonSolverList
tolPoisson               =                  1.0E-8 ! Abort criterion
abs_tol                  =                  0.0E+0 ! Lower bound for tolerance
tolCond                  =                 1.0E-23 ! Preconditioner tolerance
maxIterPoisson           =                    5000 ! Maximum iterations
poissonSolverType        =              'bicgstab' ! 'bicgstab', 'gcr', 'adi'
                                                   ! or 'hypre'
preconditioner           =                   'yes' ! 'no' or 'yes'
dtau                     =                  8.0E-1 ! Time parameter for
                                                   ! preconditioner
maxIterADI               =                      10 ! Preconditioner iterations
initialCleaning          =                       T ! Enforce initial non-
                                                   ! divergence
pressureScaling          =                       F ! Scale by P
correctMomentum          =                       T ! Correct momentum so that
                                                   ! divergence constraint is
                                                   ! fulfilled
correctDivError          =                       F ! Subtract divergence
tolcrit                  =                   'abs' ! 'abs' or 'rel'
&end

&atmosphereList
referenceQuantities      =                 'Klein' ! 'Klein', 'WKB', 'SI' or
                                                   ! 'general'
specifyReynolds          =                       F ! Use inverse Reynolds
                                                   ! number
ReInv                    =                  0.0E+0 ! Inverse Reynolds number
mu_viscous_dim           =                  0.0E+0 ! Kinematic viscosity
mu_conduct_dim           =                  0.0E+0 ! Heat conductivity
background               =            'HeldSuarez' ! 'realistic',
                                                   ! 'isothermal',
                                                   ! 'isentropic', 'const-N',
                                                   ! 'diflapse' or
                                                   ! 'HeldSuarez'
N_BruntVaisala_dim       =                  1.8E-2 ! Buoyancy frequency for
                                                   ! 'const-N'
theta0_dim               =                  3.0E+2 ! Background potential
                                                   ! temperature for
                                                   ! 'isentropic'
Temp0_dim                =                  3.0E+2 ! Background temperature
                                                   ! for 'isothermal'
press0_dim               =                  1.0E+5 ! Ground pressure
backgroundFlow_dim(1)    =                  0.0E+0 ! Initial wind
backgroundFlow_dim(2)    =                  0.0E+0 ! Initial wind
backgroundFlow_dim(3)    =                  0.0E+0 ! Initial wind
f_Coriolis_dim           =                  1.0E-4 ! Coriolis frequency
corset                   =              'constant' ! 'constant' or 'periodic'
z_tr_dim                 =                  1.2E+4 ! Tropopause height
theta_tr_dim             =                  3.0E+2 ! Potential temperature in
                                                   ! troposphere
gamma_t                  =                  0.0E+0 ! Lapse rate in troposphere
gamma_s                  =                  0.0E+0 ! Lapse rate in
                                                   ! stratosphere
tp_strato_dim            =                  2.0E+2 ! Temperature in
                                                   ! stratosphere for
                                                   ! 'HeldSuarez'
tp_srf_trp_dim           =                 3.15E+2 ! Tropical surface
                                                   ! temperature for
                                                   ! 'HeldSuarez'
tpdiffhor_tropo_dim      =                  3.0E+1 ! Temperature difference
                                                   ! between poles and
                                                   ! tropics for 'HeldSuarez'
ptdiffvert_tropo_dim     =                  2.0E+1 ! Vertical potential
                                                   ! temperature difference
                                                   ! in troposphere for
                                                   ! 'HeldSuarez'
&end

&topographyList
topography               =                       F ! Terrain-following
                                                   ! coordinates
ipolTFC                  =                       2 ! Interpolation in the
                                                   ! transformation of w
freeSlipTFC              =                       F ! Transformed free-slip
                                                   ! condition
testTFC                  =                       F ! Various TFC tests
topographyTime           =                  0.0E+0 ! Topography growth time
mountainHeight_dim       =                  5.0E+2 ! Maximum height
mountainWidth_dim        =                  1.0E+6 ! Half width
mountain_case            =                       2 ! Predefined topography
range_factor             =                  1.0E+1 ! Ratio between large and
                                                   ! small scales
&end

&boundaryList
rhoFluxCorr              =                       F ! ...
iceFluxCorr              =                       F ! ...
uFluxCorr                =                       F ! ...
vFluxCorr                =                       F ! ...
wFluxCorr                =                       F ! ...
thetaFluxCorr            =                       F ! ...
nbCellCorr               =                       1 ! ...
spongeLayer              =                       T ! General sponge layer
                                                   ! switch
sponge_uv                =                       T ! Sponge layer for
                                                   ! horizontal wind if
                                                   ! unifiedSponge = .false.
spongeHeight             =                  3.3E-1 ! Relative height of lower
                                                   ! sponge layer edge (scale
                                                   ! height for unifiedSponge
                                                   ! = .true. and spongeType
                                                   ! = 'exponential')
spongeAlphaZ_dim         =                  2.5E+2 ! Maximum relaxation rate
                                                   ! for unifiedSponge =
                                                   ! .true.
spongeAlphaZ_fac         =                  1.0E+0 ! Sponge layer factor for
                                                   ! unifiedSponge = .false.
unifiedSponge            =                       F ! Unified sponge for both
                                                   ! time schemes, applied to
                                                   ! wind and density
lateralSponge            =                       F ! Lateral sponge for
                                                   ! unifiedSponge = .true.
spongeType               =            'polynomial' ! Sponge layer profile for
                                                   ! unifiedSponge = .true.
spongeOrder              =                       1 ! Order of polynomial
                                                   ! sponge
cosmoSteps               =                       1 ! Relative strength of
                                                   ! COSMO sponge
&end

&boundaryList2
xBoundary                =              'periodic' ! Boundary conditions in x
                                                   ! ('periodic' only)
yBoundary                =              'periodic' ! Boundary conditions in y
                                                   ! ('periodic' only)
zBoundary                =            'solid_wall' ! Boundary conditions in z
                                                   ! ('periodic' or
                                                   ! 'solid_wall')
&end

&wkbList
rayTracer                =                       F ! Ray-tracer switch
&end

&tracerList
tracerSetup              =               'alpha_z' ! initial tracer
                                                   ! distribution
include_trfrc_lo         =                       T ! leading-order GW tracer
                                                   ! forcing
include_trfrc_no         =                       T ! next-order GW tracer
                                                   ! forcing
include_trfrc_mix        =                       T ! diffusive tracer mixing
&end

&iceList
inN                      =                       0 ! ...
inQ                      =                       0 ! ...
inQv                     =                       0 ! ...
nVarIce                  =                       0 ! ...
dt_ice                   =                  0.0E+0 ! ...
no_ice_source            =                       F ! ...
&end