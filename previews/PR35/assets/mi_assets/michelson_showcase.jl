using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

const cm = 1e-2
const mm = 1e-3

function reset_beamlet!(beam::GaussianBeamlet)
    BeamletOptics._drop_beams!(beam)
    r_c = first(BeamletOptics.rays(beam.chief))
    r_w = first(BeamletOptics.rays(beam.waist))
    r_d = first(BeamletOptics.rays(beam.divergence))
    BeamletOptics.intersection!(r_c, nothing)
    BeamletOptics.intersection!(r_w, nothing)
    BeamletOptics.intersection!(r_d, nothing)
    beam.chief.rays = [r_c]
    beam.waist.rays = [r_w]
    beam.divergence.rays = [r_d]
    return nothing
end

splitter_origin = [18.81cm, 23.5cm,0]

asset_dir = @__DIR__

NBK7 = DiscreteRefractiveIndex([632.8e-9], [1.51509])

## Laser
laser_assembly = MeshDummy(joinpath(asset_dir, "Laser Assembly.stl"))

# Mirror
mirror_holder = MeshDummy(joinpath(asset_dir, "Mirror Assembly.stl"))
rpm = RightAnglePrismMirror(25e-3, 25e-3)
zrotate3d!(rpm, deg2rad(45))
translate3d!(rpm, [0,33.5cm,0])

mirror_assembly = ObjectGroup([rpm, mirror_holder])

# Beamsplitter
splitter_holder = MeshDummy(joinpath(asset_dir, "Splitter Assembly.stl"))
cbs = CubeBeamsplitter(BeamletOptics.inch, NBK7)
zrotate3d!(cbs, deg2rad(-90))

splitter_assembly = ObjectGroup([cbs, splitter_holder])

# Arms
arm_holder_1 = MeshDummy(joinpath(asset_dir, "Arm Assembly 1.stl"))
arm_holder_2 = MeshDummy(joinpath(asset_dir, "Arm Assembly 2.stl"))
m1 = RoundPlanoMirror(BeamletOptics.inch, 5e-3)
zrotate3d!(m1, deg2rad(-90))
translate3d!(m1, [22cm,0,0])
m2 = RoundPlanoMirror(BeamletOptics.inch, 5e-3)
zrotate3d!(m2, deg2rad(-90))
translate3d!(m2, [12cm,0,0])

arm_1 = ObjectGroup([m1, arm_holder_1])
arm_2 = ObjectGroup([m2, arm_holder_2])

# PD
pd_holder = MeshDummy(joinpath(asset_dir, "PD Assembly.stl"))
pd = Photodetector(5e-3, 200)
translate3d!(pd, [0, -12cm, 0])

pd_assembly = ObjectGroup([pd, pd_holder])

system = System(
    [
        laser_assembly,
        mirror_assembly,
        splitter_assembly,
        arm_1,
        arm_2,
        pd_assembly
    ]
)

## setup system (dummies and optics)
translate_to3d!(mirror_assembly, [0,-10cm,0])
translate_to3d!(splitter_assembly, [18.81cm,23.5cm,0])

translate_to3d!(arm_1, splitter_origin)
translate_to3d!(arm_2, splitter_origin)
translate3d!(arm_1, [3.81cm/2, 0, 0])
translate3d!(arm_2, [0, 3.81cm/2, 0])
zrotate3d!(arm_2, deg2rad(90))

translate_to3d!(pd_assembly, splitter_origin)
translate3d!(pd_assembly, [0, -3.81cm/2, 0])

## laser fig
λ = 632.8e-9
w0 = 0.65e-3 / 2

θ_spec = 1.4e-3 / 2  # half-angle in mrad
θ_ideal = BeamletOptics.divergence_angle(λ, w0, 1)

M2 = θ_spec / θ_ideal

beam = GaussianBeamlet([0.,0,0], [0.,1,0], λ, w0; M2)

laser_fig = Figure(size=(600, 220))
laser_ax = Axis3(laser_fig[1,1], aspect=:data, elevation=0, azimuth=2π)

hidedecorations!(laser_ax)
hidespines!(laser_ax)

render!(laser_ax, beam, flen=30cm)
render!(laser_ax, laser_assembly)

save("mi_laser_assembly.png", laser_fig, px_per_unit=8)

## waist fig
zs = 0:1e-3:50cm
w, R, ψ, w0 = BeamletOptics.gauss_parameters(beam, zs)

params_fig = Figure(size=(600, 220))
waist_ax = Axis(params_fig[1,1], xlabel="z [cm]", ylabel="Beam radius [mm]")
lines!(waist_ax, zs*1e2, w * 1e3, color=RGBAf(1,0,0,1))
lines!(waist_ax, zs*1e2, -w * 1e3, color=RGBAf(1,0,0,.5), linestyle=:dashdot)
hlines!(waist_ax, 0, color=:black)
text!(waist_ax, "Optical axis")

save("mi_waist_curve.png", params_fig, px_per_unit=8)

## mirror fig
system = System([laser_assembly, mirror_assembly])
solve_system!(system, beam)

mirror_camera_view = [
    -0.383435   0.923568  -2.35922e-16  -0.105822
    -0.259272  -0.107641   0.959787      0.0303688
     0.886429   0.368016   0.280729     -0.286356
     0.0        0.0        0.0           1.0
]

fig = Figure(size=(600, 300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam)

set_view(ax, mirror_camera_view)
save("mi_corner_mirror.png", fig; px_per_unit=8, update = false)

## splitter fig
reset_beamlet!(beam)

system = System([laser_assembly, mirror_assembly, splitter_assembly])
solve_system!(system, beam)

splitter_view = [
    -0.725337   0.688394  -9.71445e-17  -0.0445975
    -0.170124  -0.179254   0.968982      0.0972088
     0.667041   0.702838   0.247132     -0.421977
     0.0        0.0        0.0           1.0
]

splitter_fig = Figure(size=(600, 400))
display(splitter_fig)
ax = LScene(splitter_fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam; flen=5cm)

set_view(ax, splitter_view)
save("mi_beamsplitter.png", splitter_fig; px_per_unit=8, update = false)

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

mi_arms_view = [
    -0.00057493  -1.0          -2.91976e-16   0.23897
    0.511649    -0.000294162   0.859195     -0.0413112
   -0.859195     0.000493977   0.511649     -0.258694
    0.0          0.0           0.0           1.0
]

fig = Figure(size=(600, 500))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam)

set_view(ax, mi_arms_view)
save("mi_arms.png", fig; px_per_unit=8, update = false)

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

pd_view = [
    -0.909654   0.415368  9.71445e-16   0.10082
    -0.188776  -0.413419  0.890757      0.100153
     0.369992   0.81028   0.45448      -0.466288
     0.0        0.0       0.0           1.0
]

fig = Figure(size=(600, 600))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)
render!(ax, beam; flen=5cm)

set_view(ax, pd_view)
save("mi_pd.png", fig; px_per_unit=8, update = false)

## fringe plot
fringes_fig = Figure(size=(600, 270))
heat1 = Axis(fringes_fig[1, 1], xlabel="x [mm]", ylabel="y [mm]", title="Before rotation", aspect=1)
heat2 = Axis(fringes_fig[1, 2], xlabel="x [mm]", ylabel="y [mm]", title="After rotation", aspect=1, yaxisposition=:right)

hidedecorations!(heat1)
hidedecorations!(heat2)

hm = heatmap!(heat1, pd.x*1e3, pd.y*1e3, intensity(pd), colormap=:viridis)

zrotate3d!(m1, 1e-3)
empty!(pd)
solve_system!(system, beam)

hm = heatmap!(heat2, pd.x*1e3, pd.y*1e3, intensity(pd), colormap=:viridis)

save("mi_fringes.png", fringes_fig, px_per_unit=4)

## phase plot
zrotate3d!(m1, -1e-3)

n = 100
Δy = λ/n

Δy * n == λ

P = zeros(n+1)

for i in eachindex(P)
    empty!(pd)
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

## intro fig
intro_camera_view = [
    -0.707107   0.707107  -5.55112e-17   0.0188074;
    -0.498749  -0.498749   0.708871      0.222826;
     0.501247   0.501247   0.705338     -0.759889;
     0.0        0.0        0.0           1.0;
]

fig = Figure(size=(600, 600))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

transparency = true
alpha = 0.1
render!(ax, beam)
render!(ax, laser_assembly; transparency, alpha)
render!(ax, mirror_holder; transparency, alpha)
render!(ax, splitter_holder; transparency, alpha)
render!(ax, arm_holder_1; transparency, alpha)
render!(ax, arm_holder_2; transparency, alpha)
render!(ax, pd_holder; transparency, alpha)

render!(ax, rpm)
render!(ax, cbs)
render!(ax, m1)
render!(ax, m2)
render!(ax, pd)

set_view(ax, intro_camera_view)
save("mi_intro_fig.png", fig; px_per_unit=8, update = false)