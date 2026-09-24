"""
```julia
set_visualization_theme!()
```

Configure Makie.jl to use a customized theme.
"""
set_visualization_theme!

function set_visualization_theme!()
    color = [:blue, :red]
    set_theme!(
        theme_latexfonts();
        fontsize = 16,
        Axis = (
            width = 200,
            height = 200,
            xgridvisible = false,
            ygridvisible = false,
            titlefont = :regular,
        ),
        Axis3 = (
            width = 200,
            height = 200,
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            titlefont = :regular,
        ),
        colormap = :seismic,
        palette = (color = color, patchcolor = color),
    )
    return
end
