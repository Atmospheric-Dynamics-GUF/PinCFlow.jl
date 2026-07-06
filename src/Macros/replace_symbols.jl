function replace_symbols end

function replace_symbols(
    x::Any,
    replacements::Vararg{Pair{Symbol, Symbol}};
    escape_others::Bool = true,
)::Any
    if x isa Expr
        for (index, arg) in enumerate(x.args)
            if arg isa Symbol
                replaced = false
                for replacement in replacements
                    if arg === replacement[1]
                        x.args[index] = replacement[2]
                        replaced = true
                    end
                end
                !replaced && escape_others && (x.args[index] = esc(arg))
            else
                x.args[index] = replace_symbols(arg, replacements...)
            end
        end
    end

    return x
end
