using GLMakie

function plot_var(semi)

    (; grid, cache) = semi
    (; var) = cache
    (; dx, dz, nx, nz, ny, zTFC) = grid

    wTFC = copy(var.w)
    for i = 1:nx
        for j = 1:ny
            for k = 0:nz
                wTFC[i, j, k] = PinCFlow_dev.vertWind(i, j, k, semi)
            end
        end
    end
    xTFC = copy(zTFC[:, 1, :])
    xTFC .= 1.0
    xTFC .= xTFC .* semi.grid.x
    @show typeof(wTFC)
    @show typeof(xTFC)
    @show typeof(zTFC)
    @show typeof(xTFC.parent)
    @show typeof(zTFC.parent)
    @show typeof(zTFC[:, 1, :].parent)
    @show typeof(wTFC[:, 1, :].parent .* semi.equations.uRef)
    fig1 = Figure()
    # ax = Axis(fig1[1, 1], title = "")
    #     hm = contourf!(
    #         ax,
    #         xTFC,
    #         zTFC[:,1,:],
    #         wTFC[:,1,:],
    #     )

    ax = Axis(fig1[1, 1], title = "")
    hm = contourf!(ax, xTFC, zTFC[:, 1, :], semi.cache.var.w[:, 1, :])
    display(fig1)

end
