using Documenter
using DocumenterMermaid
using PinCFlow

# Copy README file.
cp("README.md", "docs/src/index.md"; force = true)

# Generate documentation.
makedocs(;
    sitename = "PinCFlow Documentation",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Theory" => [
            "Continuous equations" => "theory/continous_equations.md",
            "Temporal discretization" => "theory/temporal_discretization.md",
            "Spatial discretization" => "theory/spatial_discretization.md",
        ],
        "Examples" => [
            "Mountain-wave simulation" => "examples/mountain_wave_simulation.md",
            "WKB mountain-wave simulation" => "examples/wkb_mountain_wave_simulation.md",
        ],
        "Reference" => [
            "Types" => "reference/types.md",
            "MPIOperations" => "reference/mpi_operations.md",
            "Boundaries" => "reference/boundaries.md",
            "FluxCalculator" => "reference/flux_calculator.md",
            "PoissonSolver" => "reference/poisson_solver.md",
            "Update" => "reference/update.md",
            "MSGWaM" => "reference/msgwam.md",
            "Integration" => "reference/integration.md",
            "Output" => "reference/output.md",
        ],
    ],
    pagesonly = true,
    format = Documenter.HTML(;
        repolink = "https://gitlab.dkrz.de/atmodynamics-goethe-universitaet-frankfurt/pinc.git",
        size_threshold = nothing,
        size_threshold_warn = nothing,
    ),
)
