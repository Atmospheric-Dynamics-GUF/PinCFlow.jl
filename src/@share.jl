macro share end

macro share(x::Expr)
    if x.head === :for
        return quote
            @batch $(esc(x))
        end
    else
        return quote
            @.. thread = Threaded() $(esc(x))
        end
    end
end

macro share(r::Expr, x::Expr)
    return quote
        @batch reduction = $(esc(r)) $(esc(x))
    end
end
