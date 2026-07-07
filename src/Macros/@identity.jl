macro identity end

macro identity(x::Vararg{Any})
    return x[end]
end
