using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

# Shared refractive index used for the dielectric parts below; the exact value does not
# matter for a catalog thumbnail.
n_generic = λ -> 1.5

# Catalog tiles share one camera/figure style so that they read as a set (see
# `docs/src/basics/components/components.md`): Figure(size=(400, 300)), an Axis3 with
# aspect=:data and decorations/spines hidden, at the same azimuth/elevation.

## Flat-mirror family: SquarePlanoMirror2D, SquarePlanoMirror, RectangularPlanoMirror
## (the generic `Mirror` type reuses this render too, see catalog.json).
sq2d = SquarePlanoMirror2D(15e-3)
translate3d!(sq2d, [-25e-3, 0, 0])

sq = SquarePlanoMirror(15e-3, 4e-3)

rect = RectangularPlanoMirror(20e-3, 12e-3, 4e-3)
translate3d!(rect, [25e-3, 0, 0])

flat_mirror_fig = Figure(size=(400, 300))
flat_mirror_ax = Axis3(flat_mirror_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(flat_mirror_ax)
hidespines!(flat_mirror_ax)
render!(flat_mirror_ax, sq2d)
render!(flat_mirror_ax, sq)
render!(flat_mirror_ax, rect)
autolimits!(flat_mirror_ax)
save("flat_mirror_family_showcase.png", flat_mirror_fig; px_per_unit=4, update=false)

## RightAnglePrismMirror
rapm = RightAnglePrismMirror(15e-3, 20e-3)

rapm_fig = Figure(size=(400, 300))
rapm_ax = Axis3(rapm_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(rapm_ax)
hidespines!(rapm_ax)
render!(rapm_ax, rapm)
autolimits!(rapm_ax)
save("right_angle_prism_mirror_showcase.png", rapm_fig; px_per_unit=4, update=false)

## Retroreflector
rr = Retroreflector(15e-3)

rr_fig = Figure(size=(400, 300))
rr_ax = Axis3(rr_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(rr_ax)
hidespines!(rr_ax)
render!(rr_ax, rr)
autolimits!(rr_ax)
save("retroreflector_showcase.png", rr_fig; px_per_unit=4, update=false)

## Prism / RightAnglePrism (the generic `Prism` type is built the same way)
prism = RightAnglePrism(15e-3, 20e-3, n_generic)

prism_fig = Figure(size=(400, 300))
prism_ax = Axis3(prism_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(prism_ax)
hidespines!(prism_ax)
render!(prism_ax, prism)
autolimits!(prism_ax)
save("prism_showcase.png", prism_fig; px_per_unit=4, update=false)

## ThinBeamsplitter / RoundThinBeamsplitter
tbs = RoundThinBeamsplitter(20e-3)

tbs_fig = Figure(size=(400, 300))
tbs_ax = Axis3(tbs_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(tbs_ax)
hidespines!(tbs_ax)
render!(tbs_ax, tbs)
# `RoundThinBeamsplitter` is a zero-thickness disc: the y-extent of its mesh is exactly
# zero, so `autolimits!` falls back to a large default span along that axis. Pad the
# limits manually instead so the disc reads at the same scale as the other tiles.
limits!(tbs_ax, -0.011, 0.011, -0.011, 0.011, -0.011, 0.011)
save("thin_beamsplitter_showcase.png", tbs_fig; px_per_unit=4, update=false)

## RectangularCompensatorPlate
cmp = RectangularCompensatorPlate(20e-3, 15e-3, 4e-3, n_generic)

cmp_fig = Figure(size=(400, 300))
cmp_ax = Axis3(cmp_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(cmp_ax)
hidespines!(cmp_ax)
render!(cmp_ax, cmp)
autolimits!(cmp_ax)
save("compensator_plate_showcase.png", cmp_fig; px_per_unit=4, update=false)

## PolarizationFilter
pf = PolarizationFilter(15e-3)

pf_fig = Figure(size=(400, 300))
pf_ax = Axis3(pf_fig[1, 1], aspect=:data, azimuth=0.3π, elevation=0.25π)
hidedecorations!(pf_ax)
hidespines!(pf_ax)
render!(pf_ax, pf)
# Same zero-thickness case as `RoundThinBeamsplitter` above.
limits!(pf_ax, -0.008, 0.008, -0.008, 0.008, -0.008, 0.008)
save("polarization_filter_showcase.png", pf_fig; px_per_unit=4, update=false)
