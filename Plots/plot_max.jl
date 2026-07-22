# RUN in PinCFlow/PinCFlow.jl directory:
# julia --project=examples test/plot_max.jl

using HDF5
using CairoMakie
using Statistics

# =========================
# Input / Output
# =========================
run = "3001_01"
directory = "adv"

input_file  = "/work/bb1097/b383844/PinCFlow/$directory/results/ice_mountain_wave_$(run).h5"
output_file = "../$directory/visualization/ice_mountain_max_nan_$(run).png"

# =========================
# Read data
# =========================
data = h5open(input_file, "r")

ny = size(data["q"], 2)
iy = cld(ny, 2)

tmax = size(data["q"], 4)
println("tmax = ", tmax)

t = 1:tmax

# =========================
# Fields to analyse
# =========================
fields = Dict(
    :q   => data["q"],
    :qv  => data["qv"],
    :dqv => data["iaux3"],
    :n   => data["n"],
    :w   => data["w"]
)

# =========================
# Storage
# =========================
maxvals = Dict(k => fill(NaN, tmax) for k in keys(fields))
nancnts = Dict(k => zeros(Int, tmax) for k in keys(fields))

# =========================
# Compute maxima and NaNs
# =========================
for (name, A) in fields
    for ti in t
        slice = A[:, iy, :, ti]

        finite_vals = slice[isfinite.(slice)]
        if !isempty(finite_vals)
            maxvals[name][ti] = maximum(finite_vals)
        end

        nancnts[name][ti] = count(x -> !isfinite(x), slice)
    end
end

close(data)

# =========================
# Normalize maxima to [0,1]
# =========================
for (name, vals) in maxvals
    finite_vals = vals[isfinite.(vals)]
    if !isempty(finite_vals)
        max_over_time = maximum(finite_vals)
        if max_over_time != 0.0
            vals ./= max_over_time
        end
    end
end

# =========================
# Diagnostics (optional)
# =========================
println("\nNormalized maxima ranges:")
for (name, vals) in maxvals
    finite_vals = vals[isfinite.(vals)]
    if !isempty(finite_vals)
        println(name, ": ",
            minimum(finite_vals), " – ", maximum(finite_vals))
    end
end

# =========================
# Plot
# =========================
fig = Figure(size = (1200, 900))

ax1 = Axis(fig[1, 1],
    xlabel = "t",
    ylabel = "normalized maximum",
    title  = "Zeitliche Entwicklung der Maximalwerte"
)

ax2 = Axis(fig[2, 1],
    xlabel = "t",
    ylabel = "count",
    title  = "Zeitliche Entwicklung der NaN-Werte"
)

for (name, vals) in maxvals
    lines!(ax1, t, vals;
        label  = String(name),
        linewidth = 2
    )
end

ylims!(ax1, 0.0, 1.05)
axislegend(ax1; position = :rb)

for (name, vals) in nancnts
    lines!(ax2, t, vals;
        label  = String(name),
        linewidth = 2
    )
end

axislegend(ax2; position = :rt)

save(output_file, fig)

println("\nPlot gespeichert: $output_file")
