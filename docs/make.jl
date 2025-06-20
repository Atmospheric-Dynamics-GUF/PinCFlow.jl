
using Documenter,
    PinCFlow,
    PinCFlow.Types,
    PinCFlow.PoissonSolver,
    PinCFlow.MPIOperations,
    PinCFlow.Boundaries,
    PinCFlow.Update

makedocs(;
    sitename = "My Documentation",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Namelists" => "namelists.md",
        "PoissonSolver" => "poisson.md",
    ],
)
