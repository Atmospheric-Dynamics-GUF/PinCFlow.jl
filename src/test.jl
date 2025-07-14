module submodule

struct State{A <: AbstractFloat}
    domain :: A
    grid :: A
end

function State()
    domain = 1.
    grid = 2.
    return State(domain, grid)
end

function update_variable!(state::State)
    (; domain) = state
    domain = 2.5  # Change the local variable
    return 
end

export State, update_variable!

end

include("test2.jl")
using .MoreSubs

using .submodule

st1 = State()
println(st1.domain, st1.grid)  # 
update_variable!(st1)
println(st1.domain)  # Output: 1.0 (unchanged)
