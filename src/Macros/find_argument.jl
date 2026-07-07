"""
```julia
find_argument(input::Any; argument::Any = nothing)::Any
```

If `input` is an expression, return the first occurrence of an argument to a `Val` call which itself is an argument to a function call.

If the described pattern is not found, this function returns `nothing`.

# Arguments

  - `input`: Input object.

# Keywords

  - `argument`: Argument found in `input`, to be communicated between recursive calls of this function. Note that if `find_argument` is called with `argument` set to any value other than `nothing`, it will always return that value.
"""
function find_argument end

function find_argument(input::Any; argument::Any = nothing)::Any
    if input isa Expr
        if input.head === :call && any(
            arg isa Expr &&
            arg.head === :call &&
            arg.args[1] === :Val &&
            (arg.args[2] isa Expr || arg.args[2] isa Symbol) for
            arg in input.args
        )
            for arg in input.args
                if arg isa Expr &&
                   arg.head === :call &&
                   arg.args[1] === :Val &&
                   (arg.args[2] isa Expr || arg.args[2] isa Symbol)
                    argument === nothing && (argument = arg.args[2])
                    break
                end
            end
        else
            for arg in input.args
                argument === nothing &&
                    (argument = find_argument(arg; argument))
            end
        end
    end

    return argument
end
