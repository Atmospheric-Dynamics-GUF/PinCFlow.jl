# Changelog

PinCFlow.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1) used in the Julia ecosystem. Notable changes will be documented in this file for human readability.

## Release 2.0.0

  - The auxiliary states in the wave-packet examples have been parallelized.

  - The examples now use fewer MPI processes.

  - In the documentation of the function `correct!`, the vertical indices of the buoyancy frequency were incorrect. This has been rectified.

  - Some methods of the functions `apply_blocked_layer_scheme!`, `compute_gw_integrals!`, `compute_gw_tendencies!`, and `activate_orographic_source!` did not return `nothing`, due to missing return statements. The latter have been added.

  - The mind maps in the developer guide have been updated.

  - Objects involved in PinCFlow.jl's implementation of the BiCGSTAB algorithm have been renamed to respect its correct modern spelling ("BiCGSTAB" instead of "BicGStab"). The documentation has also been adjusted accordingly.

  - A precompilation block has been added. This block precompiles the model for a one-time-step simulation in its default configuration, thus reducing the compilation time later on.

  - PinCFlow.jl's implementation of the BiCGSTAB algorithm had a convergence criterion dependent on two different averages of the residual, namely a global RMS and a global RMS of a global vertical mean. The latter prevented convergence in very specific configurations (e.g., a Boussinesq gravity-wave packet), even though it is not actually needed. It has therefore been removed from the criterion. This change has also made the code significantly more efficient.

  - Experimental features are now being tagged in the documentation.

  - There were two bugs in the computation of background fields for the following combinations of atmosphere-namelist parameters.

    ```julia
    background == LapseRates() && (troposphere_lapse_rate != 0 || stratosphere_lapse_rate != 0)
    ```

    The bugs were as follows.

     1. The computation of the mass-weighted potential temperature was actually a computation of the pressure.

     1. Two different potential temperature profiles were computed for the troposphere and stratosphere, leading to a discontinuity at the tropopause.

    Both bugs have been fixed.

## Release 1.1.1

  - A bug has been fixed where `uold` was not being assigned correctly in `update`.

## Release 1.1.0

  - Potential-temperature fluxes due to heat conduction are no longer computed in pseudo-incompressible or Boussinesq mode.

  - Two bugs that prevented MS-GWaM from being run in `SingleColumn` or `SteadyState` mode have been fixed.

  - A bug that prevented the helper function `replace_assignments` (used in the tests) from properly overwriting variables with strings has been fixed. This function now also issues a warning if an assignment wasn't found.

  - The writing of attributes to the model output file has been serialized. This fixes a bug that led to occasional HDF5 errors in parallel simulations.

  - A bug has been fixed in the averaging of orographic wavenumbers in `apply_blocked_layer_scheme!`.

  - The developer guide has been extended with information on running and updating tests, and instructions for creating new releases.

  - The behavior of the default value of the output-namelist parameter `iin` has been changed, so that it no longer results in an error but in the selection of the last record in the input file.

  - Several bugs have been fixed in the scripts for the wave-packet examples.

  - A record of the configured namelists is now included in the model output.

## Release 1.0.0

  - The documentation has been updated, corrected and improved.

  - The density reconstructions and fluxes have been removed in Boussinesq mode.

  - The sponges are now configured with functions.

      - The following sponge-namelist parameters have been removed.

          - `sponge_extent`

          - `alpharmax`

          - `betarmax`

          - `lateral_sponge`

          - `sponge_type`

          - `sponge_order`

          - `cosmo_steps`

          - `perturbation_period`

          - `perturbation_amplitude`

          - `relaxation wind`

      - The following sponge-namelist parameters have been added.

          - `lhs_sponge`

          - `rhs_sponge`

          - `relaxed_u`

          - `relaxed_v`

          - `relaxed_w`

        Each of these must be a function that takes the three spatial coordinates, the time and the time step as arguments and returns a single value. The functions `lhs_sponge` and `rhs_sponge` are used to compute the respective Rayleigh-damping coefficients, whereas the other three functions define the wind obtained through the relaxation enforced by the LHS sponge (if the parameter `relax_to_mean` is set to `false`).

  - The internal horizontal-coordinate arrays have been parallelized, making the model slightly more efficient.

  - A bug has been fixed in the construction of the `Realistic` atmospheric background.

  - A bug has been fixed in the re-dimensionalization of the GW tracer fluxes and tracer flux convergence.

  - Metadata has been added to the output, including long variable names, dimensions and labels in LaTeX format.

  - A bug has been fixed in the computation of the ray-volume-array size (this only had an impact in configurations with `wave_modes > 1`).

  - Type constraints for structure fields have been improved. Namelist parameters that previously had to be of a subtype of `AbstractFloat` now need to be of a subtype of `Real`.

  - The following atmospheric backgrounds (only available in Boussinesq mode) have been renamed.

      - `UniformBoussinesq` $\rightarrow$ `NeutralStratification`

      - `StratifiedBoussinesq` $\rightarrow$ `StableStratification`

  - The function `check_rays` no longer triggers `exit()` calls but instead raises errors. The error messages now provide more details.

  - Fixed a bug that resulted in an incorrect initialization of nonzero density fluctuations in Boussinesq mode.

  - Fixed a bug in the initialization of the mass-weighted potential temperature in compressible mode. This only had an impact when the model was initialized with nonzero density or potential-temperature fluctuations.

  - The following atmosphere-namelist parameters have been removed.

      - `initial_thetap`

    The initial potential-temperature fluctuations are now always such that $P = \overline{\rho} \overline{\theta}$.

  - The visualization function `plot_contours` of PinCFlow.jl's `CairoMakie` extension has been replaced with the function `plot_output`, with the following changes.

      - It has a simplified call signature.

      - It can also visualize ray volumes.

      - A bug has been fixed in the setting of colorbar-tick labels.

      - The space and time units are configurable.

      - Multiple variables can be visualized in one figure.

      - The colorbar labels are set automatically.

      - The background color of plots in $x$-$z$ or $y$-$z$ plane is black.

  - The following new examples have been added.

      - A 2D cold bubble

      - A 2D hot bubble

      - A 2D vortex

      - A 3D wave packet

      - A 3D WKB wave packet

    The existing examples have been modified slightly. The run and visualization scripts have been merged.

  - The default values of the domain-namelist parameters `x_size`, `y_size` and `z_size` have been changed to `1`.

## Release 0.1.0

  - First public release of PinCFlow.jl.
