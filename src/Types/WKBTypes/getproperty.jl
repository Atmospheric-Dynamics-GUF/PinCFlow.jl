@ivy @inline function getproperty(
    rays::Rays,
    name::Symbol,
)::AbstractArray{<:AbstractFloat}
    if name === :x
        getfield(rays, :data)[1, :, :, :, :]
    elseif name === :y
        getfield(rays, :data)[2, :, :, :, :]
    elseif name === :z
        getfield(rays, :data)[3, :, :, :, :]
    elseif name === :k
        getfield(rays, :data)[4, :, :, :, :]
    elseif name === :l
        getfield(rays, :data)[5, :, :, :, :]
    elseif name === :m
        getfield(rays, :data)[6, :, :, :, :]
    elseif name === :dxray
        getfield(rays, :data)[7, :, :, :, :]
    elseif name === :dyray
        getfield(rays, :data)[8, :, :, :, :]
    elseif name === :dzray
        getfield(rays, :data)[9, :, :, :, :]
    elseif name === :dkray
        getfield(rays, :data)[10, :, :, :, :]
    elseif name === :dlray
        getfield(rays, :data)[11, :, :, :, :]
    elseif name === :dmray
        getfield(rays, :data)[12, :, :, :, :]
    elseif name === :dens
        getfield(rays, :data)[13, :, :, :, :]
    else
        getfield(rays, name)
    end
end
