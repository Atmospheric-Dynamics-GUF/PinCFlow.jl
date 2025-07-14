using Documenter
using PinCFlow

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
        size_threshold_ignore = ["Types.md", "MSGWaM.md"],
    ),
)
