## intro fig
GLMakie.activate!(; ssao=true)

const intro_camera_view = [
    -0.707107   0.707107  -5.55112e-17   0.0188074;
    -0.498749  -0.498749   0.708871      0.222826;
     0.501247   0.501247   0.705338     -0.759889;
     0.0        0.0        0.0           1.0;
]

take_screenshot("mi_intro_fig.png", system, nothing; size=(600, 600), view=intro_camera_view)

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

render!(laser_ax, beam, flen=30cm)
render!(laser_ax, laser_assembly)

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

save("mi_waist_curve.png", params_fig, px_per_unit=8)

## mirror fig
system = System([laser_assembly, mirror_assembly])
solve_system!(system, beam)

const mirror_camera_view = [
    -0.383435   0.923568  -2.35922e-16  -0.105822
    -0.259272  -0.107641   0.959787      0.0303688
     0.886429   0.368016   0.280729     -0.286356
     0.0        0.0        0.0           1.0
]

take_screenshot("mi_corner_mirror.png", system, beam; view=mirror_camera_view)

## splitter fig
reset_beamlet!(beam)

system = System([laser_assembly, mirror_assembly, splitter_assembly])
solve_system!(system, beam)

splitter_fig = Figure(size=(600, 440))
splitter_ax = Axis3(splitter_fig[1,1], aspect=:data, elevation=pi/2, azimuth=2pi)

const splitter_view = [
    -0.725337   0.688394  -9.71445e-17  -0.0445975
    -0.170124  -0.179254   0.968982      0.0972088
     0.667041   0.702838   0.247132     -0.421977
     0.0        0.0        0.0           1.0
]

take_screenshot("mi_beamsplitter.png", system, beam; size=(600, 400), view=splitter_view, flen=5cm)

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

const mi_arms_view = [
    -0.00057493  -1.0          -2.91976e-16   0.23897
    0.511649    -0.000294162   0.859195     -0.0413112
   -0.859195     0.000493977   0.511649     -0.258694
    0.0          0.0           0.0           1.0
]

take_screenshot("mi_arms.png", system, beam; size=(600, 500), view=mi_arms_view, px_per_unit=8)

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

const pd_view = [
    -0.909654   0.415368  9.71445e-16   0.10082
    -0.188776  -0.413419  0.890757      0.100153
     0.369992   0.81028   0.45448      -0.466288
     0.0        0.0       0.0           1.0
]

take_screenshot("mi_pd.png", system, beam; size=(600, 600), view=pd_view, flen=5cm)

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

##
# fig = Figure()
# ax = LScene(fig[1,1])
# render!(ax, system)
# render!(ax, beam)
# BeamletOptics.render_object_normals!(ax, pd.shape)

# fig