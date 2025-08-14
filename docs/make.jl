using Documenter
using DocumenterMermaid
using PinCFlow

# Copy README file.
cp("README.md", "docs/src/index.md"; force = true)

# Copy example plots.
for file in readdir("examples/results/"; join = true)
    cp(file, "docs/src/" * file; force = true)
end

# Generate documentation.
makedocs(;
    sitename = "PinCFlow",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Mountain-wave simulation" => "examples/mountain_wave_simulation.md",
            "WKB mountain-wave simulation" => "examples/wkb_mountain_wave_simulation.md",
        ],
        "Theory" => [
            "Physics" => "theory/physics.md",
            "Numerics" => "theory/numerics.md",
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
