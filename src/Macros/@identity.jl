macro identity end

macro identity(x::Vararg{Any})
    return esc(x[end])
end
