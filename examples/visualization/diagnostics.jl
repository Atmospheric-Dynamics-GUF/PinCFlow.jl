using Pkg

Pkg.activate("examples")

using MPI
using HDF5
using CairoMakie
using Revise
using PinCFlow



function cellwidths_from_centers(x::AbstractVector)
    n = length(x)
    @assert n ≥ 2
    edges = similar(x, float(eltype(x)), n+1)
    edges[1]   = x[1] - (x[2]-x[1])/2
    edges[end] = x[end] + (x[end]-x[end-1])/2
    @inbounds for q in 2:n
        edges[q] = (x[q-1] + x[q]) / 2
    end
    return @views edges[2:end] .- edges[1:end-1]
 end

 function compute_omega_hat(
    nn::AbstractFloat,
    kpr::AbstractFloat, 
    mr::AbstractFloat,
    )::AbstractFloat

    return nn * abs(kpr) / abs(mr)

end

function compute_total_energy(
    N::Array{<: AbstractFloat, 4},
    x::AbstractVector,
    z::AbstractVector,
    k::AbstractVector,
    m::AbstractVector
)::AbstractFloat
    Δx = cellwidths_from_centers(x)
    Δz = cellwidths_from_centers(z)
    Δk = cellwidths_from_centers(k)
    Δm = cellwidths_from_centers(m[Int.(length(m)/2 +1):end])
    Δm = [reverse(Δm);Δm]

    @assert size(N,1) == length(x)
    @assert size(N,2) == length(z)
    @assert size(N,3) == length(k)
    @assert size(N,4) == length(m)
    E = 0.0
    @inbounds for ix in eachindex(x)
        wx = Δx[ix]
        for iz in eachindex(z)
            wz = Δz[iz]
            for ik in eachindex(k)
                wk = Δk[ik]
                for im in eachindex(m)
                        wm = Δm[im]
                        E +=  compute_omega_hat(nn, k[ik], m[im])* N[ix,iz,ik,im] * wx*wz*wk*wm
                end
            end
        end
    end
    return E
end

function compute_total_spectrum(
    N::Array{<: AbstractFloat, 4},
    x::AbstractVector,
    z::AbstractVector,
    k::AbstractVector,
    m::AbstractVector
)::AbstractFloat
    Δx = cellwidths_from_centers(x)
    Δz = cellwidths_from_centers(z)
    Δk = cellwidths_from_centers(k)
    Δm = cellwidths_from_centers(m[Int.(length(m)/2 +1):end])
    Δm = [reverse(Δm);Δm]

    @assert size(N,1) == length(x)
    @assert size(N,2) == length(z)
    @assert size(N,3) == length(k)
    @assert size(N,4) == length(m)
    A = 0.0
    @inbounds for ix in eachindex(x)
        wx = Δx[ix]
        for iz in eachindex(z)
            wz = Δz[iz]
            for ik in eachindex(k)
                wk = Δk[ik]
                for im in eachindex(m)
                        wm = Δm[im]
                        A +=   N[ix,iz,ik,im] * wx*wz*wk*wm
                end
            end
        end
    end
    return A
end

time_slice = 21
nn = 0.02
nt = 30
h5open("triad_wave_propagation3_new.h5") do data
        x = data["x"][:] 
        y = data["y"][:] 
        z = data["z"][1, 1, :]
        k = data["kp"][:]
        m = data["m"][:]
        t = data["t"][:]
        total_energy = zeros(nt)
        total_spectrum = zeros(nt)
        for ti in 1:nt
            wave_spectrum = data["wavespectrum"][:, 1, :, :, :, ti]

            total_energy[ti] = compute_total_energy(wave_spectrum, x, z, k, m)
            total_spectrum[ti] = compute_total_spectrum(wave_spectrum, x, z, k, m)
        end

        
        fig = Figure()
        ax1 = Axis(fig[1, 1]; xlabel="t", ylabel="Total energy")
        ax2 = Axis(fig[1, 2]; xlabel = "t", ylabel= "Total wave action")
        lines!(ax1, t[1:nt], total_energy)  # <-- attach plot to the axis in fig
        lines!(ax2, t[1:nt], total_spectrum)

        save("examples/results/total_budget_wkb_new.svg", fig)

    return
end



