"""
```julia
Nsedimentation(N::AbstractFloat)
```

Calculate the sedimentation term for the number concentration of ice crystals.

```julia
Qsedimentation(Q::AbstractFloat)
```

Calculate the sedimentation term for the mass mixing ratio of ice crystals.

```julia
Qvsource(Qv::AbstractFloat, Qv_ref::AbstractFloat)
```

Calculate the source term for water vapor based on the difference between the current water vapor mixing ratio and a reference profile (initial profile).
"""

function Nsedimentation(N::AbstractFloat, tau_sink::AbstractFloat)

    return - 0.5 * N / tau_sink
end

function Qsedimentation(Q::AbstractFloat, tau_sink::AbstractFloat)

    return - 1.0 * Q / tau_sink
end

function Qvsource(Qv::AbstractFloat, Qv_ref::AbstractFloat, tau_relax::AbstractFloat)

    return 1.0 * (Qv_ref - Qv) / tau_relax
end

# complex sink term: qsink = - 1.0 / tau * q[i, j, k]^(5. / 3.) * n[i, j, k]^(-2. / 3.)
# complex sink term: nsink = - 0.5 / tau * q[i, j, k]^(2. / 3.) * n[i, j, k]^(1. / 3.)