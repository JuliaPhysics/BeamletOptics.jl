# Wedge light trap

A wedge light trap folds a laser beam many times between two surfaces at a small angle, a mirror and an
edge filter whose cut-off wavelength depends on the angle of incidence. It was proposed to multiply the
phase modulation that a moving mirror imposes on a beam [Ewers:2019](@cite). The beam enters through the
filter at a high angle of incidence, at which the filter transmits. Every round trip in the wedge reduces
the angle of incidence by twice the wedge angle, until the filter reflects: the beam is trapped. It runs
towards the apex of the wedge, turns at normal incidence and runs back, until the filter lets it out
again near the point where it entered.

## The edge filter

An interference filter shifts its edge to shorter wavelengths at oblique incidence. For a filter of the
effective index ``n_\mathrm{eff}`` and the edge ``\lambda_0`` at normal incidence,

```math
\lambda_\mathrm{edge}(\theta) = \lambda_0 \sqrt{1 - \frac{\sin^2\theta}{n_\mathrm{eff}^2}},
```

with ``\theta`` the angle of incidence in air. A longpass filter transmits a laser of the wavelength
``\lambda`` if ``\lambda > \lambda_\mathrm{edge}(\theta)``. A [`SpectralCoating`](@ref) takes this as a
function of the wavelength and the angle of incidence in air; the coating sees the same Snell invariant
from the glass of its substrate and from the air of the wedge.

```@example wolit
using GLMakie, BeamletOptics
GLMakie.activate!(; ssao = true) # hide
const BMO = BeamletOptics
const mm = 1e-3

λ = 632.8e-9                  # HeNe laser
λ_cut, n_eff = 650e-9, 1.85   # edge at normal incidence, effective index of the filter
edge(θ) = λ_cut * sqrt(1 - sin(θ)^2 / n_eff^2)

longpass = SpectralCoating((λ, θ) -> λ > edge(θ) ? 0.0 : 1.0; splitting = BMO.Dominant())
nothing # hide
```

The filter transmits the HeNe line above an angle of incidence of

```@example wolit
rad2deg(asin(n_eff * sqrt(1 - (λ / λ_cut)^2)))
```

degrees.

## The wedge

Two glass plates form the wedge, 1 mm apart at the origin and tilted by ±0.6°. The lower plate carries an
anti-reflection coating on its outer face and the filter on its inner face, the upper plate a mirror on
its inner face. The faces of a plate are addressed by their outward normals in the local frame of the
plate, see [Coatings](@ref).

```@example wolit
α = deg2rad(1.2)   # wedge angle
glass = ConstantRefractiveIndex(1.5)

filter_plate = RectangularCompensatorPlate(40mm, 10mm, 1mm, glass)
surface_model!(filter_plate, FlatFace([0, -1, 0]), AntireflectionCoating())
surface_model!(filter_plate, FlatFace([0, 1, 0]), longpass)
translate3d!(filter_plate, [0, -1mm, 0])
zrotate3d!(filter_plate, α / 2)

mirror_plate = RectangularCompensatorPlate(40mm, 10mm, 1mm, glass)
surface_model!(mirror_plate, FlatFace([0, -1, 0]), MirrorCoating())
translate3d!(mirror_plate, [0, 1mm, 0])
zrotate3d!(mirror_plate, -α / 2)

system = System([filter_plate, mirror_plate])
nothing # hide
```

## Trapping the beam

The beam enters at 26.5°, just above the edge of the filter.

```@example wolit
θ_in = deg2rad(26.5)
beam = Beam([-12mm, -4mm, 0], [sin(θ_in), cos(θ_in), 0], λ)
solve_system!(system, beam; r_max = 400)

fig = Figure(size = (1000, 360))
ax = LScene(fig[1, 1]; show_axis = false)
render!(ax, system)
render!(ax, beam; color = :red, flen = 4mm)
set_view(ax, [-8mm, 0.5mm, 14mm], [-8mm, 0.5mm, 0], [0, 1, 0]) # hide
set_orthographic(ax) # hide
fig
```

The angle of incidence at the mirror drops by ``2\alpha = 2.4°`` per round trip, turns at normal
incidence and rises again, until the beam leaves through the filter:

```@example wolit
in_wedge(r) = BMO.refractive_index(r) == 1 && BMO.direction(r)[2] > 0
aoi = [rad2deg(acos(BMO.direction(r)[2])) for r in rays(beam) if in_wedge(r)]

fig = Figure(size = (600, 340))
ax = Axis(fig[1, 1]; xlabel = "round trip", ylabel = "angle of incidence [°]")
scatterlines!(ax, aoi)
hlines!(ax, rad2deg(asin(n_eff * sqrt(1 - (λ / λ_cut)^2))); color = :gray, linestyle = :dash)
fig
```

## The acceptance window

Only beams whose angle of incidence falls below the edge within one round trip are trapped: the trap
accepts a window of ``2\alpha`` above the edge.

```@example wolit
θs = deg2rad.(24.5:0.1:28.5)
round_trips = map(θs) do θ
    b = Beam([-12mm, -4mm, 0], [sin(θ), cos(θ), 0], λ)
    solve_system!(system, b; r_max = 1000)
    count(in_wedge, rays(b))
end

fig = Figure(size = (600, 340))
ax = Axis(fig[1, 1]; xlabel = "angle of incidence [°]", ylabel = "round trips in the wedge")
scatterlines!(ax, rad2deg.(θs), round_trips)
fig
```

Below the edge the filter reflects the beam before it enters the wedge.
