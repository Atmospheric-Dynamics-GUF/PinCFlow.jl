l2 = (
    buoyancyproduction = 0.0021437851f0,
    kek = 6.7082014f0,
    kh = 6.7082014f0,
    km = 6.7082014f0,
    n2 = 0.010105865f0,
    p = 5973.465f0,
    pip = 4.4238324f-5,
    rhobar = 17.034365f0,
    rhop = 0.0013121589f0,
    shearproduction = 1.355889f-9,
    t = 3600.0f0,
    thetabar = 13598.672f0,
    tke = 0.002236067f0,
    us = 447.23724f0,
    vs = 0.32979524f0,
    wts = 0.544529f0,
    x = 18165.902f0,
    y = 18165.902f0,
    z = 364725.28f0,
    ztilde = 392459.12f0,
)
linf = (
    buoyancyproduction = 6.8040375f-5,
    kek = 0.21217465f0,
    kh = 0.21217465f0,
    km = 0.21217465f0,
    n2 = 0.00031957557f0,
    p = 321.18332f0,
    pip = 1.0270111f-5,
    rhobar = 1.0363004f0,
    rhop = 0.00045842316f0,
    shearproduction = 4.3861675f-10,
    t = 3600.0f0,
    thetabar = 556.882f0,
    tke = 5.002105f-5,
    us = 10.127179f0,
    vs = 0.047865972f0,
    wts = 0.13445848f0,
    x = 9000.0f0,
    y = 9000.0f0,
    z = 19001.666f0,
    ztilde = 20000.0f0,
)
reference = (l2, linf)

@testset "Mountain wave" begin
    test_example(
        joinpath(scripts_directory, "mountain_wave.jl"),
        reference,
        :x_size => 10,
        :y_size => 10,
        :z_size => 10,
        :output => OutputNamelist(; prepare_restart = true);
        update_references,
    )
end
