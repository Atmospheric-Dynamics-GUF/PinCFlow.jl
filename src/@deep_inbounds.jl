macro deep_inbounds end

macro deep_inbounds(x::Expr)
    return quote
        $(deep_inbounds(esc(x)))
    end
end
