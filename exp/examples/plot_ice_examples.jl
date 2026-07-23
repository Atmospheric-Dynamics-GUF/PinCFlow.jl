
#import the modules
using HDF5
#using Statistics
#using PyPlot: contourf, colorbar, xlabel, ylabel, color
using Plots
using LaTeXStrings
pyplot()

file_path = "ice_mountain_2D_lapserates.h5"


function plot_nucleation(file_path::String; mode::Symbol=:both, t_idx::Int=7)
    # Read grid dimensions and time
    x = h5read(file_path, "x") .* 0.001        # convert to km
    z_raw = h5read(file_path, "z")
    z = vec(z_raw[1, 1, :]) .* 0.001          # convert to km
    t = h5read(file_path, "t")
    current_time = t[t_idx]

    # Helper function to read transposed 2D slice (x, z)
    read_slice(var_name) = h5read(file_path, var_name)[:, 1, :, t_idx]'

    # load wind components
    u = read_slice("u")
    w = read_slice("w")

    if mode == :both
        # --- 3 x 2 Layout (6 plots) ---
        n_hom = read_slice("n_hom")
        n_het = read_slice("n_het")
        q_hom = read_slice("q_hom")
        q_het = read_slice("q_het")

        p1 = contourf(x, z, n_hom, title=L"$n_{hom}$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$n_{hom}$ [1/kg]", color=:blues, levels=20, lw=0)
        p2 = contourf(x, z, n_het, title=L"$n_{het}$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$n_{het}$ [1/kg]", color=:blues, levels=20, lw=0)
        p3 = contourf(x, z, q_hom, title=L"$q_{hom}$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$q_{hom}$ [kg/kg]", color=:blues, levels=20, lw=0)
        p4 = contourf(x, z, q_het, title=L"$q_{het}$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$q_{het}$ [kg/kg]", color=:blues, levels=20, lw=0)
        p5 = contourf(x, z, u,     title=L"$u$ (t=%$current_time s)",     xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$u$ [m/s]",     color=:viridis, levels=20, lw=0)
        p6 = contourf(x, z, w,     title=L"$w$ (t=%$current_time s)",     xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$w$ [m/s]",     color=:vik,     levels=20, lw=0)

        return plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(900, 950), margin=4Plots.mm, aspect_ratio=:auto)

    elseif mode == :hom || mode == :het
        # --- 2 x 2 Layout (4 plots) ---
        is_hom = (mode == :hom)
        n_var = is_hom ? "n_hom" : "n_het"
        q_var = is_hom ? "q_hom" : "q_het"
        lbl   = is_hom ? "hom"   : "het"

        n_data = read_slice(n_var)
        q_data = read_slice(q_var)

        p1 = contourf(x, z, n_data, title=LaTeXString("\$n_{$(lbl)}\$ (t=$(current_time) s)"), xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=LaTeXString("\$n_{$(lbl)}\$ [1/kg]"), color=:blues, levels=20, lw=0)
        p2 = contourf(x, z, q_data, title=LaTeXString("\$q_{$(lbl)}\$ (t=$(current_time) s)"), xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=LaTeXString("\$q_{$(lbl)}\$ [kg/kg]"), color=:blues, levels=20, lw=0)
        p3 = contourf(x, z, u,      title=L"$u$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$u$ [m/s]", color=:viridis, levels=20, lw=0)
        p4 = contourf(x, z, w,      title=L"$w$ (t=%$current_time s)", xlabel=L"$x$ [km]", ylabel=L"$z$ [km]", colorbar_title=L"$w$ [m/s]", color=:vik, levels=20, lw=0)

        return plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 700), margin=4Plots.mm, aspect_ratio=:auto)
    else
        error("Invalid mode! Use :both, :hom, or :het")
    end
end

# configure files
files_to_process = [
    (path="ice_mountain_2D_test_bothnucl.h5", mode=:both),
    (path="ice_mountain_2D_test_homnucl.h5", mode=:hom),
    (path="ice_mountain_2D_test_hetnucl.h5", mode=:het)
]

# Generate and display figures
for item in files_to_process
    fig = plot_nucleation(item.path, mode=item.mode, t_idx=7)
    display(fig)
end