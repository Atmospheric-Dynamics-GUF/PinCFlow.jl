using Documenter
using PinCFlow

# Copy README file.
cp("README.md", "docs/src/index.md"; force = true)

# Generate documentation.
makedocs(;
    sitename = "PinCFlow Documentation",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Reference" => [
            "Types.md",
            "MPIOperations.md",
            "Boundaries.md",
            "FluxCalculator.md",
            "PoissonSolver.md",
            "Update.md",
            "MSGWaM.md",
            "Integration.md",
            "Output.md",
        ],
    ],
    pagesonly = true,
    format = Documenter.HTML(;
        repolink = "git@gitlab.dkrz.de:atmodynamics-goethe-universitaet-frankfurt/pinc.git",
        size_threshold = nothing,
        size_threshold_warn = nothing,
    ),
)
