function ivy end

function ivy(x::Any)::Any
    if isa(x, Expr)
        head = string(x.head)
        args = string.(x.args)

        if head in ("=", ".=")
            x.args[2] = ivy(x.args[2])
        elseif endswith(head, r"[^!=]=")
            assignment = replace(head, r"[^.=]" => "")
            operation = head[1:(end - 1)]
            x = Meta.parse(
                "$(args[1]) $assignment $(args[1]) $operation ($(args[2]))",
            )
            x.args[2] = ivy(x.args[2])
        elseif head == "ref"
            for index in eachindex(x.args)
                x.args[index] = ivy(x.args[index])
            end
            x = quote
                @inbounds @views $x
            end
        else
            for index in eachindex(x.args)
                x.args[index] = ivy(x.args[index])
            end
        end
    end

    return x
end
