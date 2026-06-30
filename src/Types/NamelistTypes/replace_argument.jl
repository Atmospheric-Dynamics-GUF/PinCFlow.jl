"""
```julia
replace_argument(input::Any, argument::Any, literal::Any)::Any
```

Return a deep copy of `input` with all occurrences of `argument` as argument of a `Val` call replaced by `literal` (if `input` is an expression).

# Arguments

  - `input`: Input object.

  - `argument`: Argument in `Val` calls.

  - `literal`: Replacement for `argument` in `Val` calls.
"""
function replace_argument end

function replace_argument(input::Any, argument::Any, literal::Any)::Any
    output = deepcopy(input)

    if output isa Expr
        if output.head === :call &&
           output.args[1] === :Val &&
           output.args[2] === argument
            output.args[2] = literal
        else
            for (index, arg) in enumerate(output.args)
                output.args[index] = replace_argument(arg, argument, literal)
            end
        end
    end

    return output
end
