abstract type AbstractSpongeLayer end
abstract type UnifiedSponge <: AbstractSpongeLayer end

struct ExponentialSponge <: UnifiedSponge end
struct SinusodialSponge{RealT<:Real} <: UnifiedSponge
    ## TODO: add types
    lateral::Any
    height::Any
    alphaZ_dim::Any
    alpha::AbstractArray
    relax_to_mean::RealT
    relaxation_period::RealT
    relaxation_amplitude::RealT
    kr_sp_tfc::AbstractArray
    kr_sp_w_tfc::AbstractArray
end

function SinusodialSponge(pars::Parameters)
    nx, ny, nz = pars.domain.sizex, pars.domain.sizey, pars.domain.sizez
    alpha = OffsetArray(zeros(Float64, 0:(nx+1), 0:(ny+1), 0:(nz+1)))
    kr_sp_tfc = OffsetArray(zeros(nx + 2, ny + 2, nz + 2), 0:(nx+1), 0:(ny+1), 0:(nz+1))
    kr_sp_w_tfc = copy(kr_sp_tfc)
    SinusodialSponge(true, 0.5, 0.01, alpha, 1.0, 0.0, 0.0, kr_sp_tfc, kr_sp_w_tfc)
end

struct CosmoSponge <: UnifiedSponge end
struct PolynomialSponge <: UnifiedSponge end
struct NonUnifiedSponge <: AbstractSpongeLayer end

initialize_sponge!(model) = initialize_sponge_ini!(model, model.atmosphere.sponge)

function setup_sponge(sponge_type::String, pars::Parameters)
    if sponge_type == "sinusodial"
        return SinusodialSponge(pars)
    else
        error("invalid sponge type $sponge_type")
    end
end

function initialize_sponge_ini!(model, sponge::SinusodialSponge)
    @trixi_timeit timer() "Initialize sponge" begin
    #! format: noindent
    (; nx, ny, nz) = model.domain
    (; x, y, ztfc, lx, ly, lz) = model.grid
    (; tref) = model.constants
    sponge = model.atmosphere.sponge
    alphaZ = sponge.alphaZ_dim * tref

    dz = sponge.height * (lz[1] - lz[0])
    z = lz[1] - dz

    if sponge.lateral == true
        dx = 0.5 * sponge.height * (lx[1] - lx[0])
        dy = 0.5 * sponge.height * (ly[1] - ly[0])
        x0 = lx[0] + dx
        x1 = lx[1] - dx
        y0 = ly[0] + dy
        y1 = lx[1] - dy

        if nx > 1 && ny > 1
            alphaZ = alphaZ / 3.0
            alphaX = alphaZ
            alphaY = alphaZ
        elseif nx > 1
            alphaZ = alphaZ / 2.0
            alphaX = alphaZ
            alphaY = 0.0
        elseif ny > 1
            alphaZ = alphaZ / 2.0
            alphaX = 0.0
            alphaY = alphaZ
        end
    end

    for k in 1:nz
        for j in 0:(ny+1)
            for i in 0:(nx+1)
                height = ztfc[i, j, k]

                if nz > 1
                    if height >= z # this is zSponge
                        sponge.alpha[i, j, k] += alphaZ *
                                                 sin(0.5 * pi * (height - z) / dz)^2
                    end
                end

                # todo: express this via the type system
                if sponge.lateral == true
                    if nx > 1
                        if x[i] <= x0
                            sponge.alpha[i, j, k] += alphaX *
                                                     sin(0.5 * pi * (x0 - x[i]) / dx)^2
                        elseif x[i] >= x1
                            sponge.alpha[i, j, k] += alphaX *
                                                     sin(0.5 * pi * (x[i] - x1) / dx)^2
                        end
                    end
                end
                if ny > 1
                    if y[j] <= y0
                        sponge.alpha[i, j, k] += alphaY *
                                                 sin(0.5 * pi * (y0 - y[j]) / dy)^2
                    elseif y[i] >= y1
                        sponge.alpha[i, j, k] += alphaY *
                                                 sin(0.5 * pi * (y[j] - y1) / dy)^2
                    end
                end
            end
        end
        # TODO - Use broadcasting and @views
        sponge.alpha[:, :, 0] = sponge.alpha[:, :, 1]
        sponge.alpha[:, :, nz+1] = sponge.alpha[:, :, nz]
    end
    end # timer
end
