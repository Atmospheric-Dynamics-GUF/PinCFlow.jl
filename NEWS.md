# Changelog

PinCFlow.jl follows the interpretation of
[semantic versioning (semver)](https://julialang.github.io/Pkg.jl/dev/compatibility/#Version-specifier-format-1) used in the Julia ecosystem. Notable changes will be documented in this file for human readability.

## Release 5.0.0

New features:

  - An experimental elastic-mode-selection algorithm, adapted from [Banerjee (2026)](https://doi.org/10.5281/zenodo.20582010), has been implemented. It can be used to reduce the number of ray volumes launched by the orographic source ([#194](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/194)).

  - Output fields for WKB integrals have been added ([#236](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/236))

Improvements:

  - The composite type `Preconditioner` has been removed from `Types.PoissonTypes`. The corresponding field of `Types.PoissonTypes.Poisson` has also been removed ([#220](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/220)). Since the previous release, the arrays needed by the preconditioner are in `state.variables.auxiliaries`.

  - The WKB-namelist parameter `multiplication_factor` has been replaced with the parameters `k_bins`, `l_bins`, and `m_bins`. These specify the exact number of bins in $k$, $l$, and $m$ in the ray-volume merging algorithm. Merging is now triggered when the ray-volume count in a grid cell exceeds `k_bins * l_bins * m_bins` ([#221](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/221)).

  - The `@ivy` macro now turns off all bounds checks when applied to any expression (including, e.g., function definitions) ([#226](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/226)).

  - The Makie.jl extension has been improved in several ways ([#239](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/239)).

  - The examples have been modified and made cheaper. The example plots are no longer stored in `examples/results/` but generated in the documentation workflow ([#240](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/240), [#254](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/254), [#261](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/261)).

  - The macros `@ivy` and `@dispatch` have been moved to their own module (`Macros`), and `@dispatch` has been made more robust - it can no longer modify expressions or strings by accident ([#241](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/241)).

  - The reduction of exceptions in MPI simulations now works across nodes. The function `reduce_exceptions` has been moved to the `Integration` module ([#243](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/243)).

  - The function `ensemble` has been removed. Instead, ensemble simulations can now be run by simply passing multiple namelists objects to `integrate` ([#243](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/243)).

  - A few expensive runtime allocations in MS-GWaM have been removed ([#250](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/250))

  - The field `rhs` of `Types.PoissonTypes.Poisson` has been renamed `lhs` ([#252](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/252))

  - The calculation of turbulent velocities in MS-GWaM has been improved ([#253](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/253)).

  - Refraction has been switched off for ray volumes that were launched in the previous time step. The function `compute_vertical_averages` in `MSGWaM.Raysources` has been replaced with the new function `compute_orographic_flow` in `MSGWaM.BlockedLayer`. The blocked-layer scheme has been improved ([#244](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/244)).

  - Several errors in the documentation have been fixed, and it has been improved in a few places ([#225](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/225), [#241](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/241), [#258](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/258), [#261](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/261)).

Bug fixes:

  - MS-GWaM can no longer be initialized with multiple ray-volumes in single-grid-cell dimensions, i.e., an error is thrown for `(x_size == 1 && nrx > 1) || (y_size == 1 && nry > 1)` ([#221](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/221)).

  - Several bugs related to the TKE scheme have been fixed ([#228](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/228), [#234](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/234), [#238](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/238), [#261](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/261)).

  - A bug in the ray-volume launch algorithm of the orographic source has been fixed ([#229](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/229)).

  - A minor bug in the calculation of WKB integrals has been fixed ([#246](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/246)).

  - A bug in the calculation of heat conduction has been fixed ([#248](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/248))

  - Tracer-related bugs in the I/O and `smooth_gw_tendencies!` have been fixed ([#258](https://github.com/Atmospheric-Dynamics-GUF/PinCFlow.jl/pull/258)).

## Release 4.0.0

  - PinCFlow.jl has been extended with a turbulence parameterization. Using the new prognostic variable turbulent kinetic energy (TKE), vertical diffusion can now be applied to momentum, potential temperature, and tracers. Configuration is provided through the `TurbulenceNamelist`. The turbulence parameterization is enabled by default.

  - An experimental feature has been added that damps the wave-action density of ray volumes through wave-induced turbulence. As this feature has not yet been thoroughly validated, it is disabled by default. It can be enabled by setting the WKB-namelist parameter `turbulent_damping = true`.

  - This release also includes two bug fixes:

      - The averaging using the box filter has been adjusted to respect the irregular grid.

      - A minor correction to the mean-flow interpolation. This bug had a negligible impact over uneven terrain and no impact over flat terrain.

  - These's been a minor correction to the documentation of 'Grid'.

## Release 3.0.0

Breaking changes, performance improvements, and compilation acceleration:

  - PinCFlow.jl now requires the Julia version 1.10.11 or higher.

  - To avoid a naming conflict with CairoMakie's `Box`, the WKB-namelist parameter `filter_type` now takes `:BoxFilter` and `:ShapiroFilter` instead of `:Box` and `:Shapiro`, respectively.

  - The user can now define a vertical grid stretch function via the grid-namelist parameter `vertical_grid_stretching`.

  - The initial potential-temperature fluctuations can be defined by the user via the atmosphere-namelist parameter `initial_thetap` when `buoyancy_initialization == true`.

  - The following changes have been made to the tracer transport:

      - With future extensions in mind, the gravity-wave-tracer fluxes and tendencies have been restructured by removing `TracerWKBImpact` and `TracerForcing` and adding `TracerWKBIntegrals` and `TracerWKBTendencies`, to match the structures of `WKBIntegrals` and `WKBTendencies`.

      - When setting the tracer-namelist parameter `apply_sponge_to_tracer == true`, the tracer fields are relaxed to the user-defined `relaxed_chi`.

      - The tracer is now initialized via the tracer-namelist parameter `initial_chi`.

  - Namelist parameters of singleton types have been replaced with parameters of the type `Symbol` and namelist parameters that are functions are now stored in `FunctionWrapper`s to avoid massive recompilation due to a changed concrete `State` type.

  - Project activation has been removed from all scripts, with the project now specified via the `--project` flag.

  - The precompilation block now runs cheap versions of the examples.

  - The communication of ray volumes has been made more efficient.

  - The number of MPI processes used in the example shell scripts has been increased.

  - The function `ensemble` for conducting ensemble runs has been added.

  - The example scripts have been replaced by functions that are exported by the `PinCFlow` module.

  - The example shell scripts for the Goethe cluster have been removed.

  - A setup script has been added to configure the `examples` project's backends on the Levante cluster and precompile, which can be triggered by sourcing the file `examples/levante/setup.sh`.

Bug fixes and further minor improvements:

  - There was a bug in the calculation of the stratospheric mass-weighted potential temperature for the atmospheric-namelist parameter `background == :Realistic`, which has been fixed.

  - There were a few dependency issues with Revise.jl, CairoMakie.jl, HDF5.jl and MPI.jl which have been fixed.

  - The handling of the propagation of ray volumes has been improved to allow for scenarios with steep topography and/or high vertical resolution.

  - A bug has been fixed in the blocked-layer scheme, and it has been moved to its own submodule. The orographic source has undergone minor restructuring.

  - A few checks have been added to ensure proper parallelization initialization.

  - The `create_output` function now creates the output directory if it doesn't exist yet.

  - The `integrate` function now prints the approximate peak memory usage across all MPI processes.

  - The function `reduce_exceptions` has been added to handle exceptions across MPI processes.

Documentation corrections:

  - The list of publications has been updated.

  - There were minor errors in the documentation of gravity-wave fluxes.

  - Minor aesthetic adjustments have been made to the documentation.

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

    The initial potential-temperature fluctuations are now always such that $P = \bar{\rho} \bar{\theta}$.

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
