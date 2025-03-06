## intro fig
intro_fig = Figure(size=(600, 380))
intro_ax = Axis3(intro_fig[1,1], aspect=:data, elevation=0.4, azimuth=8.7)

hidedecorations!(intro_ax)
hidespines!(intro_ax)

render_system!(intro_ax, system)

save("mi_intro_fig.png", intro_fig, px_per_unit=4)

## laser fig
λ = 632.8e-9
w0 = 0.65e-3 / 2

θ_ideal = BeamletOptics.divergence_angle(λ, w0, 1) * 1e3

M2 = 1.4 / θ_ideal

beam = GaussianBeamlet([0.,0,0], [0.,1,0], λ, w0; M2)

laser_fig = Figure(size=(600, 220))
laser_ax = Axis3(laser_fig[1,1], aspect=:data, elevation=0, azimuth=2π)

hidedecorations!(laser_ax)
hidespines!(laser_ax)

render_beam!(laser_ax, beam, flen=30cm)
render_object!(laser_ax, laser_assembly)

save("mi_laser_assembly.png", laser_fig, px_per_unit=4)

## waist fig
zs = 0:1e-3:50cm
w, R, ψ, w0 = BeamletOptics.gauss_parameters(beam, zs)

params_fig = Figure(size=(600, 220))
waist_ax = Axis(params_fig[1,1], xlabel="z [cm]", ylabel="Beam radius [mm]")
lines!(waist_ax, zs*1e2, w * 1e3, color=RGBAf(1,0,0,1))
lines!(waist_ax, zs*1e2, -w * 1e3, color=RGBAf(1,0,0,.5), linestyle=:dashdot)
hlines!(waist_ax, 0, color=:black)
text!(waist_ax, "Optical axis")

save("mi_waist_curve.png", params_fig, px_per_unit=4)

## mirror fig
system = System([laser_assembly, mirror_assembly])
solve_system!(system, beam)

mirror_fig = Figure(size=(600, 300))
mirror_ax = Axis3(mirror_fig[1,1], aspect=:data, elevation=0.4, azimuth=2π)

hidedecorations!(mirror_ax)
hidespines!(mirror_ax)

render_system!(mirror_ax, system)
render_object!(mirror_ax, rpm)
render_beam!(mirror_ax, beam, flen=10cm)
render_object!(mirror_ax, laser_assembly)

save("mi_corner_mirror.png", mirror_fig, px_per_unit=4)

## splitter fig
reset_beamlet!(beam)

system = System([laser_assembly, mirror_assembly, splitter_assembly])
solve_system!(system, beam)

splitter_fig = Figure(size=(600, 440))
splitter_ax = Axis3(splitter_fig[1,1], aspect=:data, elevation=pi/2, azimuth=2pi)

hidedecorations!(splitter_ax)
hidespines!(splitter_ax)

render_system!(splitter_ax, system)
render_object!(splitter_ax, cbs)
render_beam!(splitter_ax, beam, flen=10cm)
render_object!(splitter_ax, rpm)
render_object!(splitter_ax, laser_assembly)

save("mi_beamsplitter.png", splitter_fig, px_per_unit=4)

## arms fig
reset_beamlet!(beam)

system = System([
    laser_assembly,
    mirror_assembly,
    splitter_assembly,
    arm_1,
    arm_2
])
solve_system!(system, beam)

arm_fig = Figure(size=(600, 390))
arm_ax = Axis3(arm_fig[1,1], aspect=:data, elevation=0.45, azimuth=5.3)

hidedecorations!(arm_ax)
hidespines!(arm_ax)

render_system!(arm_ax, system)
render_beam!(arm_ax, beam, flen=20cm)

save("mi_arms.png", arm_fig, px_per_unit=4)

## system figure
reset_beamlet!(beam)

system = System([
    laser_assembly,
    mirror_assembly,
    splitter_assembly,
    arm_1,
    arm_2,
    pd_assembly
])

solve_system!(system, beam)

pd_fig = Figure(size=(600, 550))
pd_ax = Axis3(pd_fig[1,1], aspect=:data, elevation=pi/2, azimuth=2pi)

hidedecorations!(pd_ax)
hidespines!(pd_ax)

render_system!(pd_ax, system)
render_object!(pd_ax, cbs)
render_beam!(pd_ax, beam)
render_object!(pd_ax, laser_assembly)
render_object!(pd_ax, rpm)
render_object!(pd_ax, m1)
render_object!(pd_ax, m2)

save("mi_pd.png", pd_fig, px_per_unit=4)

## fringe plot
fringes_fig = Figure(size=(600, 270))
heat1 = Axis(fringes_fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", title="Before rotation", aspect=1)
heat2 = Axis(fringes_fig[1, 2], xlabel="x [mm]", ylabel="y [mm]", title="After rotation", aspect=1, yaxisposition=:right)

hidedecorations!(heat1)
hidedecorations!(heat2)

hm = heatmap!(heat1, pd.x*1e3, pd.y*1e3, BeamletOptics.intensity(pd), colormap=:viridis)

zrotate3d!(m1, 1e-3)
BeamletOptics.reset_detector!(pd)
solve_system!(system, beam)

hm = heatmap!(heat2, pd.x*1e3, pd.y*1e3, BeamletOptics.intensity(pd), colormap=:viridis)

save("mi_fringes.png", fringes_fig, px_per_unit=4)

## phase plot
zrotate3d!(m1, -1e-3)

n = 100
Δy = λ/n

Δy * n == λ

P = zeros(n+1)

for i in eachindex(P)
    BeamletOptics.reset_detector!(pd)
    solve_system!(system, beam)
    P[i] = BeamletOptics.optical_power(pd)
    # translate by Δy
    translate3d!(m2, [0, Δy, 0])
end

ys = LinRange(0, n*Δy, n+1)

power_fig = Figure(size=(600, 250))
power_ax = Axis(power_fig[1, 1], xlabel="Δy [nm]", ylabel="P [mW]",)

lines!(power_ax, ys*1e9, P*1e3, color=:red)
vlines!(power_ax, λ*1e9, color=:red, linestyle=:dashdot)

ylims!(power_ax, 0, 1)

save("mi_powerplot.png", power_fig, px_per_unit=4)

##
# fig = Figure()
# ax = LScene(fig[1,1])
# render_system!(ax, system)
# render_beam!(ax, beam)
# BeamletOptics.render_object_normals!(ax, pd.shape)

# fig