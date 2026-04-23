function getproperty(rays::Rays, name::Symbol)::AbstractArray{<:AbstractFloat}
    @ivy if name == :x
        return getfield(rays, :data)[1, :, :, :, :]
    elseif name == :y
        return getfield(rays, :data)[2, :, :, :, :]
    elseif name == :z
        return getfield(rays, :data)[3, :, :, :, :]
    elseif name == :k
        return getfield(rays, :data)[4, :, :, :, :]
    elseif name == :l
        return getfield(rays, :data)[5, :, :, :, :]
    elseif name == :m
        return getfield(rays, :data)[6, :, :, :, :]
    elseif name == :dxray
        return getfield(rays, :data)[7, :, :, :, :]
    elseif name == :dyray
        return getfield(rays, :data)[8, :, :, :, :]
    elseif name == :dzray
        return getfield(rays, :data)[9, :, :, :, :]
    elseif name == :dkray
        return getfield(rays, :data)[10, :, :, :, :]
    elseif name == :dlray
        return getfield(rays, :data)[11, :, :, :, :]
    elseif name == :dmray
        return getfield(rays, :data)[12, :, :, :, :]
    elseif name == :dens
        return getfield(rays, :data)[13, :, :, :, :]
    else
        return getfield(rays, name)
    end
end
