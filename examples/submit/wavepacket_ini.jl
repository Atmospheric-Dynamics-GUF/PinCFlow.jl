# isothermal atmosphere 

p0 = 1.0E+5 # ground pressure 
t0 = 300.0 # ground temperature
g = 9.81
r = 287
href = r * t0 / g
kappa = 2.0 / 7.0 # R/c_p

a0 = 0.5

sigmax = 0.0
sigmay = 0.0
sigmaz = 5000.0

lambdax = 30.0E+3
lambday = 0.0
lambdaz = 3000.0

x0 = 0.0
y0 = 0.0
z0 = 10.E+3

n2 = g^2.0 * kappa / t0 / r

f = 0.0
f2 = f^2.0
branch = -1.0

function pbar(z)
    return p0 * exp(-z / href)
end

function thetabar(z)
    return t0 * (p0 / pbar(z))^kappa
end

function rhobar(z)
    return pbar(z) / (r * t0)
end

function envelope(x, y, z)
    deltax = x - x0
    deltay = y - y0
    deltaz = z - z0

    if sigmax == 0.0
        envel = 1.0
    elseif abs(deltax) < sigmax
        envel = 0.5 * (1.0 + cos(deltax * pi / sigmax))
    else
        envel = 0.0
    end

    if sigmay == 0.0
        envel = 1.0 * envel
    elseif abs(deltay) < sigmay
        envel = envel * 0.5 * (1.0 + cos(deltay * pi / sigmay))
    else
        envel = 0.0
    end

    envel = envel * exp(-deltaz^2.0 / (2.0 * sigmaz^2.0))

    return envel
end

if lambdax == 0.0
    kk = 0.0
else
    kk = 2 * pi / lambdax
end

if lambday == 0.0
    ll = 0.0
else
    ll = 2 * pi / lambday
end

mm = 2 * pi / lambdaz

kh2 = kk^2.0 + ll^2.0
omega = branch * sqrt((n2 * kh2 + f2 * mm^2) / (kh2 + mm^2.0))
omega2 = omega^2.0

function bhat(x, y, z)
    return a0 * n2 / mm * envelope(x, y, z)
end

function phi(x, y, z)
    return kk * x + ll * y + mm * z
end

function bprime(x, y, z)::AbstractFloat
    return real(bhat(x, y, z) * exp(1im * phi(x, y, z)))
end

function uprime(x, y, z)::AbstractFloat
    uhat =
        1im / (mm * n2) * (omega2 - n2) / (omega2 - f2) *
        (kk * omega + 1im * ll * f) *
        bhat(x, y, z)

    return real(uhat * exp(1im * phi(x, y, z)))
end

function vprime(x, y, z)::AbstractFloat
    vhat =
        1im / (mm * n2) * (omega2 - n2) / (omega2 - f2) *
        (ll * omega - 1im * kk * f) *
        bhat(x, y, z)

    return real(vhat * exp(1im * phi(x, y, z)))
end

function wprime(x, y, z)::AbstractFloat
    what = 1im * omega / n2 * bhat(x, y, z)

    return real(what * exp(1im * phi(x, y, z)))
end

function piprime(x, y, z)::AbstractFloat
    pihat =
        kappa / r / thetabar(z) * 1im / mm * (omega2 - n2) / n2 * bhat(x, y, z)

    return real(pihat * exp(1im * phi(x, y, z)))
end

function rhoprime(x, y, z)::AbstractFloat
    return rhobar(z) / (1 + bprime(x, y, z) / g) - rhobar(z)
end

function thetaprime(x, y, z)::AbstractFloat
    return bprime(x, y, z) * thetabar(z) / g
end