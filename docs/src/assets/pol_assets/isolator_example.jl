using GLMakie, CairoMakie, BeamletOptics

const BMO = BeamletOptics

GLMakie.activate!(; ssao=true)

include(joinpath(@__DIR__, "..", "render_utils.jl"))

mm = 1e-3

wp = HalfWaveplate(12mm)
cbs = PolarizingCubeBeamsplitter(20mm, \lambda -> 1.5)

translate3d!(wp, [0, 10mm, 0])
translate3d!(cbs, [0, 30mm, 0])

beam = PolarizedRay([0,0,0], [0,1,0], 532e-9, [1,0,0])
laser = Beam(beam)

system = StaticSystem([wp, cbs])

solve_system!(system, laser)

# System render
view = [0.794 0.607 0.0 -0.04; 0.0 0.0 1.0 0.0; 0.607 -0.794 0.0 -0.15; 0 0 0 1]

take_screenshot("pol_system.png", system, laser; size=(600,400), view)

# Power vs rotation
pd_t = Photodetector(5mm, 25)
translate3d!(pd_t, [0, 50mm, 0])

pd_r = Photodetector(5mm, 25)
translate3d!(pd_r, [40mm, 30mm, 0])
zrotate3d!(pd_r, deg2rad(90))

sim = StaticSystem([wp, cbs, pd_t, pd_r])

angles = range(0, stop=deg2rad(90), length=50)
P_t = zeros(length(angles))
P_r = zeros(length(angles))

for (i,a) in enumerate(angles)
    reset_rotation3d!(wp)
    zrotate3d!(wp, a)
    empty!(pd_t); empty!(pd_r)
    solve_system!(sim, laser)
    P_t[i] = BMO.optical_power(pd_t)
    P_r[i] = BMO.optical_power(pd_r)
end

fig = Figure(size=(600,250))
ax = Axis(fig[1,1], xlabel="θ [deg]", ylabel="Relative power")
lines!(ax, rad2deg.(angles), P_t ./ maximum(P_t), color=:blue, label="transmitted")
lines!(ax, rad2deg.(angles), P_r ./ maximum(P_t), color=:red, label="reflected")
axislegend(ax)

save("pol_power.png", fig, px_per_unit=4)

