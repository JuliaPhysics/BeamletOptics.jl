using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const mm = 1e-3

l = 1e-3
R1 = 100e-3
R2 = Inf
d = 25.4e-3
n = 1.5
λ = 1e-6

cs = UniformDiscSource([0, -10mm, 0], [0, 1, 0], 15e-3, λ)
lens = SphericalLens(R1, R2, l, d, x -> n)

psfd = PSFDetector(10e-3)
translate3d!(psfd, [0, 200mm + 0.13mm, 0])

sys = System([lens, psfd])

solve_system!(sys, cs)

x, z, I_num = intensity(psfd; n=500, crop_factor=10)

fringes_fig = Figure()
srf = Axis3(fringes_fig[1, 1], xlabel="x [µm]", ylabel="z [µm]", elevation=0.9, azimuth=2.372)
hm = surface!(srf, x*1e6, z*1e6, I_num, colormap=:viridis)

hidezdecorations!(srf)

##
k = -0.675
d = 75.0e-3
l = 15e-3
radius = 76.68e-3
A = [0*(1e3)^1, 2.7709219e-8*(1e3)^3, 6.418186e-13*(1e3)^5, -1.5724014e-17*(1e3)^7, -2.7768768e-21*(1e3)^9, -2.590162e-25*(1e3)^11]
AL75150 = Lens(
    EvenAsphericalSurface(radius, d, k, A),
    l,
    n -> 1.5006520430
)

xrotate3d!(AL75150, deg2rad(-0.5))

pd = PSFDetector(15e-3)

translate3d!(pd, [0, 158.1779e-3, 0.0])
system = System([AL75150, pd])

ps = UniformDiscSource([0, -0.1, 0], [0,1,0], 0.8*d, 1550e-9)

solve_system!(system, ps)

asph_fig = Figure(size=(800, 500))

xs, ys, _intensity = intensity(pd, n=801, x_min=-25e-6, x_max=25e-6, z_min=-25e-6, z_max=25e-6, z0_shift=-50e-6)

ax = Axis(asph_fig[1,1], xlabel="x [µm]", ylabel="z [µm]", aspect=1)

hm = heatmap!(xs*1e6, ys*1e6, _intensity)
