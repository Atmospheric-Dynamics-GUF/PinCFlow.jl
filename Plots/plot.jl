# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/plot.jl

# plot given variable from h5 files created by ice_mountain_2D.jl

using HDF5
using PyPlot
using Statistics

# Base font-size (points). Change this to scale all plot text.
FS = 12

# tau values for the sink, period values for u relaxation
tau_vals = [0.0, 0.00000003, 0.0, 0.00000003]
period_vals = [1800.0, 1800.0, 18000.0, 18000.0]

# variables to plot
varnames = ["u", "w", "q"]

for varname in varnames
    figure(figsize=(12, 4 * length(tau_vals)))

    for (i, tau) in enumerate(tau_vals)
        period = period_vals[i]
        filename = "../ice_mountain_wave_tau_$(tau)_period_$(period).h5"
        data = h5open(filename, "r")

        x = data["x"][:] .* 0.001
        z = data["z"][1, 1, :] .* 0.001

        ny = size(data[varname], 2)
        iy = cld(ny, 2)

        tidx = size(data[varname], 4)
        fld = data[varname][:, iy, :, tidx]

        subplot(length(tau_vals), 1, i)
        contour = pcolormesh(x, z, fld; cmap="viridis")
        colorbar(contour)
        title("tau=$(tau), period=$(period)")
        ylabel("z (km)")
        if i == length(tau_vals)
            xlabel("x (km)")
        end

    suptitle("Field: $(varname)", fontsize=1.5*FS)
    savefig("test_output/mountain_wave_$(varname).png", dpi=150)
    close()
    end
end


