# examples/submit/wave_packet.jl

using Pkg

Pkg.activate("examples")

using Revise
using PinCFlow

npx = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
npy = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
npz = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

x_size = 40
y_size = 40
z_size = 40

lx = 20000.0
ly = 20000.0
lz = 20000.0

rx = 0.0
ry = 0.25
rz = 0.25

x0 = 0.0
y0 = 0.0
z0 = 3 * lz / 4

a0 = 0.1

k = 16 * pi / lx
l = 16 * pi / ly
m = -16 * pi / lz

background = Realistic()
coriolis_frequency = 0.0001

atmosphere = AtmosphereNamelist(; background, coriolis_frequency)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz)
auxiliary_state = State(Namelists(; atmosphere, domain))
(; g, kappa, rsp) = auxiliary_state.constants

function ijk(x, y, z)
    i = argmin(abs.(x .- auxiliary_state.grid.x))
    j = argmin(abs.(y .- auxiliary_state.grid.y))
    k = argmin(abs.(z .- auxiliary_state.grid.zc[i, j, :]))

    return CartesianIndex(i, j, k)
end

function n2(x, y, z)
    return auxiliary_state.atmosphere.n2[ijk(x, y, z)]
end

function rhobar(x, y, z)
    return auxiliary_state.atmosphere.rhobar[ijk(x, y, z)]
end

function thetabar(x, y, z)
    return auxiliary_state.atmosphere.thetabar[ijk(x, y, z)]
end

function envelope(x, y, z)
    r =
        sqrt(
            (rx * k * (x - x0))^2 +
            (ry * l * (y - y0))^2 +
            (rz * m * (z - z0))^2,
        ) / pi
    if r <= 1
        return (1 + cos(pi * r)) / 2
    else
        return 0.0
    end
end

function phi(x, y, z)
    return k * x + l * y + m * z
end

function omega(x, y, z)
    return -sqrt(
        (n2(x, y, z) * (k^2 + l^2) + coriolis_frequency^2 * m^2) /
        (k^2 + l^2 + m^2),
    )
end

function bhat(x, y, z)
    return a0 * n2(x, y, z) / m * envelope(x, y, z)
end

function uhat(x, y, z)
    return 1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (k * omega(x, y, z) + 1im * l * coriolis_frequency) *
           bhat(x, y, z)
end

function vhat(x, y, z)
    return 1im / m / n2(x, y, z) * (omega(x, y, z)^2 - n2(x, y, z)) /
           (omega(x, y, z)^2 - coriolis_frequency^2) *
           (l * omega(x, y, z) - 1im * k * coriolis_frequency) *
           bhat(x, y, z)
end

function what(x, y, z)
    return 1im * omega(x, y, z) / n2(x, y, z) * bhat(x, y, z)
end

function pihat(x, y, z)
    return kappa / rsp / thetabar(x, y, z) * 1im / m *
           (omega(x, y, z)^2 - n2(x, y, z)) / n2(x, y, z) * bhat(x, y, z)
end

atmosphere = AtmosphereNamelist(;
    background,
    coriolis_frequency,
    initial_rhop = (x, y, z) ->
        -rhobar(x, y, z) / g * real(bhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_u = (x, y, z) -> real(uhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_v = (x, y, z) -> real(vhat(x, y, z) * exp(1im * phi(x, y, z))),
    initial_w = (x, y, z) -> real(what(x, y, z) * exp(1im * phi(x, y, z))),
    initial_pip = (x, y, z) ->
        real(pihat(x, y, z) * exp(1im * phi(x, y, z))),
)
domain = DomainNamelist(; x_size, y_size, z_size, lx, ly, lz, npx, npy, npz)
output = OutputNamelist(;
    output_variables = (:u, :v, :w),
    output_file = "wave_packet.h5",
    tmax = 60.0,
    output_interval = 60.0,
)

integrate(Namelists(; atmosphere, domain, output))
