"""
```julia
replace_assignments(
    code::AbstractString,
    assignments::Vararg{Pair{Symbol, <:Any}},
)::AbstractString
```

Return `code` with replaced `assignments`.

# Arguments

  - `code`: String of Julia code.

  - `assignments`: Replacements for assignments in the code.
"""
function replace_assignments end

function replace_assignments(
    code::AbstractString,
    assignments::Vararg{Pair{Symbol, <:Any}},
)::AbstractString
    for assignment in assignments
        (name, value) = assignment
        typeof(value) <: AbstractString && (value = "\"$value\"")
        range = findfirst(Regex("$name *= *"), code)
        if range !== nothing
            (start, stop) = extrema(range)
            suffix = ""
            stop += 1
            try
                stop = Meta.parse(code, stop)[2]
                suffix = "\n"
            catch
                while !(code[stop] in (',', ')'))
                    stop = Meta.parse(code, stop; greedy = false)[2]
                end
            end
            stop -= 1
            code = replace(code, code[start:stop] => "$name = $value" * suffix)
        else
            println("Warning: No assignment of \"$name\" was found!")
            println("")
        end
    end

    return code
end
