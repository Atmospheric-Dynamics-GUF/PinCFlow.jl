macro last end

macro last(x::Vararg{Any})
    return esc(x[end])
end
