using GLMakie, BeamletOptics
using LinearAlgebra: dot, normalize

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

# The custom OpenRAMAN components (ReflectiveGrating, LongpassDichroicMirror, glasses)
# live in the `OpenRaman` module. Guarded and loaded relative to `@__MODULE__` so this
# script is self-sufficient wherever it runs: standalone (this becomes `Main.OpenRaman`),
# or reached from the tutorial's `@setup raman` block via `conditional_include` (this
# becomes a module scoped to `DocUtils`, distinct from -- but just as usable as -- the
# tutorial's own copy).
isdefined(@__MODULE__, :OpenRaman) || include(joinpath(@__DIR__, "OpenRaman", "OpenRaman.jl"))
using .OpenRaman

const mm = 1e-3
const nm = 1e-9

asset_dir = @__DIR__

# Excitation and Raman-shifted wavelengths used throughout
const λ_pump = 532nm
const λ_raman = (561nm, 588nm, 633nm)

#=
Optomechanics

The baseplate and cover are imported as `MeshDummy`s: they are rendered but ignored by the
solver, so beam clearance can be inspected in the same scene as the optical path. Every
position below is a CAD coordinate relative to the baseplate origin.
=#
baseplate = MeshDummy(joinpath(asset_dir, "Baseplate.stl"))
cover = MeshDummy(joinpath(asset_dir, "Cover.stl"))

#=
Cuvette holder

The holder carries a cemented doublet that focuses the pump into the sample and a
cylindrical lens on the entrance side. The three parts are bundled in an `ObjectGroup` so
they move as one, but only the doublet is handed to the solver -- see the tutorial for why
the cylindrical lens is held out of the traced set.
=#
cuvette_mesh = MeshDummy(joinpath(asset_dir, "Cuvette Holder.stl"))

AC127_019 = SphericalDoubletLens(12.9mm, -11mm, -59.3mm, 4.5mm, 1.5mm, 12.7mm, N_BAF10, N_SF6HT)

LK1085L1 = Lens(CylindricalSurface(-10.3mm, 15mm, 17mm), 2mm, N_BK7)
yrotate3d!(LK1085L1, deg2rad(-90))
zrotate3d!(LK1085L1, deg2rad(180))
translate_to3d!(LK1085L1, [0, 15.857mm + thickness(LK1085L1), 0])

cuvette_holder = ObjectGroup([AC127_019, LK1085L1, cuvette_mesh])
translate_to3d!(cuvette_holder, [-20.075mm, 154.648mm, 21.779mm])

#=
Excitation path: laser -> fold mirror -> dichroic -> cuvette
=#
PF10G01 = RoundPlanoMirror(25.4mm, 6mm)
translate_to3d!(PF10G01, [29.925mm, 110.721mm, 21.778mm])
zrotate3d!(PF10G01, deg2rad(-45))

DMLP550 = LongpassDichroicMirror(25.4mm, 3mm, UV_Fused_Silica, 550nm)
translate_to3d!(DMLP550, position(PF10G01) .+ [-50mm, 0, 0])
zrotate3d!(DMLP550, deg2rad(-45 + 180))

#=
Collection path: longpass filter -> window -> collimator -> grating -> focusing lens -> sensor
=#

# longpass filter
FELH0550 = Prism(BMO.PlanoSurfaceSDF(3.5mm, 21.1mm), UV_Fused_Silica)
translate_to3d!(FELH0550, [-20.941mm, 62.662mm, 21.778mm])
zrotate3d!(FELH0550, deg2rad(180 + 15))

# precision window
WG41050 = Prism(BMO.PlanoSurfaceSDF(5mm, 25.4mm), UV_Fused_Silica)
translate_to3d!(WG41050, [-20.645mm, 44.162mm, 21.778mm])
zrotate3d!(WG41050, deg2rad(180 + 22))

AC127_019_col = SphericalDoubletLens(12.9mm, -11mm, -59.3mm, 4.5mm, 1.5mm, 12.7mm, N_BAF10, N_SF6HT)
translate_to3d!(AC127_019_col, [-19.934mm, 24.162mm, 21.779mm])
zrotate3d!(AC127_019_col, deg2rad(180))

AC254_050_1 = SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, 25.4mm, N_BAF10, N_SF10)
translate_to3d!(AC254_050_1, position(AC127_019_col) .+ [0, -76.297mm, 0])

const groove_density = 1.2e6        # 1200 lines/mm
const grating_order = 1
const grating_alpha = 47.5          # mount angle in deg

GR25_1205 = RectangularReflectiveGrating(25mm, 25mm, 6mm, groove_density, grating_order)
translate_to3d!(GR25_1205, position(AC254_050_1) .+ [0, -50mm, 0])
zrotate3d!(GR25_1205, deg2rad(180 - grating_alpha))

const f_focus = 50mm                # nominal focal length of the AC254-050 focusing lens

AC254_050_2 = SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, 25.4mm, N_BAF10, N_SF10)
translate_to3d!(AC254_050_2, position(GR25_1205) + BMO.orientation(GR25_1205)[:,2] * -45mm)
zrotate3d!(AC254_050_2, deg2rad(-grating_alpha))

pd_size = 20mm
pd = Detector(pd_size)
translate_to3d!(pd, position(AC254_050_2) + BMO.orientation(AC254_050_2)[:,2] * 54.5mm)
zrotate3d!(pd, deg2rad(-grating_alpha))

#=
The system. Note that only `AC127_019` of the cuvette group is handed over: the
cylindrical lens is rendered but not traced, and the holder mesh is a `MeshDummy`.
=#
optical_system = StaticSystem([
    AC127_019, PF10G01, DMLP550, FELH0550, WG41050,
    AC127_019_col, AC254_050_1, GR25_1205, AC254_050_2, pd
])

#=
Excitation beam: a coherent 532 nm source, modelled as a single `GaussianBeamlet`
=#
laser = GaussianBeamlet(
    [29.925mm, 79.721mm, 21.778mm],
    [0, 1, 0],
    λ_pump,
    3.5mm/2
)

solve_system!(optical_system, laser)

# Where does the pump waist sit inside the cuvette?
zs = LinRange(100mm, 175mm, 1000)
w, R, ψ, ~ = gauss_parameters(laser, zs)

w_min, i_min = findmin(w)
z_focus = zs[i_min]
focal_point = point_on_beam(laser, z_focus)[1]

@info "Pump waist inside the cuvette" w_min z_focus

#=
Raman return: incoherent emission from the illuminated sample volume, modelled as
`PointSource` ray bundles -- one per Raman line. All three use the same spread angle and
ray count so that the only difference between them is the wavelength.
=#
raman_origin = position(AC127_019) .+ [0, 21.62mm, 0]
raman_angle = deg2rad(6)
raman_rays = 60
raman_rings = 3

#=
The three fans are otherwise identical -- same origin, same direction, same spread angle
and ray count -- and `PointSource` derives its sampling basis from `dir` deterministically.
Without intervention all three would place their rays at exactly the same azimuths and draw
on top of each other in the render. Passing a `basis` rotated about the fan axis spins each
one against the others, which changes nothing physically (the sampled cone is the same) but
makes the three colours distinguishable wherever they share a path.
=#
raman_axis = [0, -1, 0]
raman_basis(i) = BMO.rotate3d(normalize(raman_axis), (i - 1) * deg2rad(10)) * [1.0, 0, 0]

raman_sources = map(enumerate(λ_raman)) do (i, λ)
    PointSource(raman_origin, raman_axis, raman_angle, λ;
        num_rays=raman_rays, num_rings=raman_rings, basis=raman_basis(i))
end

spots = map(raman_sources) do src
    empty!(pd)
    solve_system!(optical_system, src)
    spot_diagram(pd)
end

raman_561, raman_588, raman_633 = raman_sources
sd_561, sd_588, sd_633 = spots

#=
Spectral calibration sweep: for a range of wavelengths, record where the ray bundle lands
on the sensor. A coarser bundle is enough here, since only the centroid is used. The sweep
deliberately starts below the dichroic cut-on, where nothing reaches the sensor at all.
=#
centroid(pts) = sum(pts) / length(pts)

λ_sweep = LinRange(540nm, 675nm, 46)

x_sweep = map(λ_sweep) do λ
    src = PointSource(raman_origin, raman_axis, raman_angle, λ; num_rays=40, num_rings=2)
    empty!(pd)
    solve_system!(optical_system, src)
    # below the DMLP550 cut-on the light never leaves the excitation arm
    isnothing(BMO.hits(pd)) && return NaN
    return centroid(spot_diagram(pd))[1]
end

#=
Analytical comparison: first-order linear dispersion of the spectrograph. The grating
equation sin(β) = mλ/d - sin(α) gives the diffraction angle, and the focusing lens converts
angle into position with dx/dλ = m·f / (d·cos β).

The remaining sign is pure bookkeeping: a `Detector` reports its hits in a *left-handed*
(x, z) frame, so that the spot diagram reads the right way round when viewed from the
incoming beam. Its local x-axis is therefore `-orientation(pd)[:,1]`.
=#
sinβ(λ) = grating_order * λ * groove_density - sind(grating_alpha)

pd_local_x = -BMO.orientation(pd)[:,1]
x_sign = sign(dot(pd_local_x, BMO.orientation(GR25_1205)[:,1]))

λ_ref = λ_raman[2]
β_ref = asin(sinβ(λ_ref))
dxdλ = x_sign * grating_order * f_focus * groove_density / cos(β_ref)

i_ref = argmin(abs.(λ_sweep .- λ_ref))
x_ref = x_sweep[i_ref]
x_analytic = x_ref .+ dxdλ .* (λ_sweep .- λ_sweep[i_ref])

@info "Linear dispersion" dxdλ*1e-9*1e6 β_ref*180/π

#=
=================================== Figures ===================================
=#

# Camera pose recorded interactively for the full-instrument views
system_view = (
    [0.14269811319617923, 0.17484439806107677, 0.294120884914244],
    [0.004158535711631309, 0.04787331363451436, 0.010519756043704757],
    [-0.6145438743625663, -0.5632275109364454, 0.5523681726573785]
)

render_all!(ax) = begin
    render!(ax, optical_system)
    render!(ax, LK1085L1)
    render!(ax, cuvette_mesh; transparency=true, alpha=0.1)
    render!(ax, laser; color=RGBAf(0,1,0,1), flen=26mm)
    # Seen from the full-instrument distance the three fans pile up into a grey blob on
    # the collection arm, before the grating disperses them. Show only every eighth ray
    # and lift the opacity so the surviving lines stay recognisably blue, green and red.
    # The grating close-up renders the fans separately at a denser sampling.
    render!(ax, raman_561; color=RGBAf(0,0,1,0.5), flen=50mm, render_every=4)
    render!(ax, raman_588; color=RGBAf(0,1,0,0.5), flen=50mm, render_every=3)
    render!(ax, raman_633; color=RGBAf(1,0,0,0.5), flen=50mm, render_every=4)
end

## intro figure -- the whole instrument
intro_fig = Figure(size=(700, 500))
display(intro_fig)
ax = LScene(intro_fig[1,1])
hide_axis(ax)

render!(ax, baseplate; transparency=true, alpha=0.1)
render!(ax, cover; transparency=true, alpha=0.05)
render_all!(ax)

set_view(ax, system_view...)
save("or_intro_fig.png", intro_fig; px_per_unit=8, update=false)

## coordinate frame figure -- where the CAD numbers come from
frame_fig = Figure(size=(700, 450))
display(frame_fig)
ax = LScene(frame_fig[1,1])
hide_axis(ax)

render!(ax, baseplate; transparency=true, alpha=0.02)
render_lcs!(ax, baseplate; scale=6, show_labels=true)
render!(ax, optical_system)
render!(ax, laser; color=RGBAf(0,1,0,1), flen=26mm)

system_view2 = [
     -0.574701   0.818363  7.7083e-9  -0.0377639
 -0.506236  -0.355508  0.785709    0.0282435
  0.642996   0.451548  0.618596   -0.286405
  0.0        0.0       0.0         1.0
]

set_view(ax, system_view2)
save("or_coordinates.png", frame_fig; px_per_unit=8, update=false)

## excitation path -- laser, fold mirror, dichroic, cuvette
exc_fig = Figure(size=(700, 450))
display(exc_fig)
ax = LScene(exc_fig[1,1])
hide_axis(ax)

render!(ax, baseplate; transparency=true, alpha=0.05)
render!(ax, cuvette_mesh; transparency=true, alpha=0.1)
render!(ax, PF10G01)
render!(ax, DMLP550)
render!(ax, AC127_019)
render!(ax, LK1085L1)
render!(ax, laser; color=RGBAf(0,1,0,1), flen=26mm)
scatter!(ax, focal_point; color=:red, markersize=18)

look_at!(ax, [0.004, 0.118, 0.022], [0.105, -0.115, 0.085])
save("or_excitation.png", exc_fig; px_per_unit=8, update=false)

## pump waist along the excitation path
waist_fig = Figure(size=(600, 250))
waist_ax = Axis(waist_fig[1,1], xlabel="Path length z [mm]", ylabel="Beam radius [mm]")
lines!(waist_ax, zs*1e3, w*1e3, color=RGBAf(1,0,0,1))
lines!(waist_ax, zs*1e3, -w*1e3, color=RGBAf(1,0,0,0.5), linestyle=:dashdot)
hlines!(waist_ax, 0, color=:black)
vlines!(waist_ax, z_focus*1e3, color=:blue, linestyle=:dashdot)
text!(waist_ax, z_focus*1e3, 0.9*maximum(w)*1e3;
    text=" waist $(round(w_min*1e6, digits=1)) µm", align=(:left, :center), fontsize=12)

save("or_waist.png", waist_fig, px_per_unit=4)

## dichroic close-up -- the pump reflects, the Raman light passes through
dic_fig = Figure(size=(700, 450))
display(dic_fig)
ax = LScene(dic_fig[1,1])
hide_axis(ax)

render!(ax, DMLP550)
render!(ax, PF10G01)
render!(ax, FELH0550)
render!(ax, laser; color=RGBAf(0,1,0,1), flen=26mm)
render!(ax, raman_633; color=RGBAf(1,0,0,0.4), flen=50mm, render_every=2)

look_at!(ax, position(DMLP550), [0.02, -0.02, 0.10])
save("or_dichroic.png", dic_fig; px_per_unit=8, update=false)

## grating close-up -- one bundle in, three colours out
gr_fig = Figure(size=(700, 450))
display(gr_fig)
ax = LScene(gr_fig[1,1])
hide_axis(ax)

render!(ax, GR25_1205)
render!(ax, AC254_050_1)
render!(ax, AC254_050_2)
render!(ax, pd)
render!(ax, raman_561; color=RGBAf(0,0,1,0.5), flen=50mm, render_every=2)
render!(ax, raman_588; color=RGBAf(0,1,0,0.5), flen=50mm, render_every=2)
render!(ax, raman_633; color=RGBAf(1,0,0,0.5), flen=50mm, render_every=2)

look_at!(ax, position(GR25_1205) .+ [0.01, 0.03, 0.0], [0.0, -0.02, 0.13])
save("or_grating.png", gr_fig; px_per_unit=8, update=false)

## spot diagram -- three Raman lines on the sensor, full sensor and per-spot zoom
spot_colors = (:blue, :green, :red)
spot_labels = ("561 nm", "588 nm", "633 nm")

spot_fig = Figure(size=(800, 400))
spot_ax = Axis(spot_fig[1,1], xlabel="x [mm]", ylabel="y [mm]", aspect=1, title="Full sensor")
zoom_ax = Axis(spot_fig[1,2], xlabel="x - x₀ [µm]", ylabel="y - y₀ [µm]", aspect=1,
    title="Centred on each spot", yaxisposition=:right)

for (sd, c, l) in zip(spots, spot_colors, spot_labels)
    scatter!(spot_ax, first.(sd)*1e3, last.(sd)*1e3; color=c, markersize=5, label=l)
    # re-centre each cluster on its own centroid to compare the geometric blur
    c0 = centroid(sd)
    scatter!(zoom_ax, (first.(sd) .- c0[1])*1e6, (last.(sd) .- c0[2])*1e6;
        color=c, markersize=5, label=l)
end

xlims!(spot_ax, -pd_size/2*1e3, pd_size/2*1e3)
ylims!(spot_ax, -pd_size/2*1e3, pd_size/2*1e3)
xlims!(zoom_ax, -60, 60)
ylims!(zoom_ax, -60, 60)
axislegend(zoom_ax, position=:rt, framevisible=false, labelsize=10)

save("or_spots.png", spot_fig, px_per_unit=4)

## spectral calibration -- simulated vs. analytical linear dispersion
cal_fig = Figure(size=(600, 300))
cal_ax = Axis(cal_fig[1,1], xlabel="Wavelength [nm]", ylabel="Sensor position x [mm]")

vlines!(cal_ax, DMLP550.cuton*1e9; color=:black, linestyle=:dot)
text!(cal_ax, DMLP550.cuton*1e9, minimum(filter(!isnan, x_sweep))*1e3;
    text=" dichroic cut-on", align=(:left, :bottom), fontsize=10)
lines!(cal_ax, λ_sweep*1e9, x_analytic*1e3; color=:blue, linestyle=:dashdot, linewidth=2,
    label="Linear dispersion\n($(round(dxdλ*1e-9*1e6, digits=1)) µm/nm)")
scatter!(cal_ax, λ_sweep*1e9, x_sweep*1e3; color=:red, markersize=6, label="Simulated")
axislegend(cal_ax, position=:rt, framevisible=false, labelsize=10)

save("or_calibration.png", cal_fig, px_per_unit=4)
