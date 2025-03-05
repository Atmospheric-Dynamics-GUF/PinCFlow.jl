using GLMakie

function plot_var(semi)
    (; grid, cache) = semi
    (; var) = cache
    (; dx, dz, nx, nz, ny) = grid

    wTFC = copy(var.w)
    for i in 1:nx
        for j in 1:ny
            for k in 0:nz
                wTFC[i, j, k] = PinCFlow_dev.vertWind(i, j, k, semi)
            end
        end
    end

    fig1 = Figure()
    ax = Axis(fig1[1, 1], title = "")
    hm = contourf!(ax,
                   grid.x[1:(grid.nx)],
                   grid.z[0:(grid.nz)] .+ dz / 2,
                   wTFC[1:(grid.nx), 1, 0:(grid.nz)] .* semi.equations.uRef)
    Colorbar(fig1[1, 2], hm)
    display(fig1)

    @show maximum(abs.(wTFC .* semi.equations.uRef))
end
