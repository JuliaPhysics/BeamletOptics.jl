using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

distance = 20e-2
factor = 1.2
RoC = distance/2 * factor
m1 = ConcaveSphericalMirror(RoC, 5e-3, 2BeamletOptics.inch)
m2 = ConcaveSphericalMirror(RoC, 5e-3, 2BeamletOptics.inch)

zrotate3d!(m1, deg2rad(180))
translate3d!(m2, [0, distance, 0])

system = StaticSystem([m1, m2])

fig = Figure(size=(600,240))
dr = 0.03
y1 = -0.02
y2 = 0.22
aspect = (1, (y2-y1)/(2dr), 1)
limits = (-dr, dr, y1, y2, -dr, dr)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0.0, elevation=1e-3)

hidedecorations!(ax)
hidespines!(ax)

beam = Beam(Ray([0, distance/2, 7e-3], [0.17, 1, 0]))

solve_system!(system, beam, r_max=100)

render!(ax, beam, flen=0.1)
render!(ax, m1)
render!(ax, m2)

save("concave_mirror_showcase.png", fig; px_per_unit=8)