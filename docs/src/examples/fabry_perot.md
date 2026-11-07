# Fabry-Pérot and Fizeau interferometers

Two partial mirrors facing each other form a cavity. The light that enters it is reflected back and
forth, and at every reflection a part of it leaves. All these parts interfere: in a Fabry-Pérot etalon
with parallel mirrors as a sharp resonance in the transmitted spectrum, between two mirrors at a small
angle as the fringes of a Fizeau interferometer. In BMO each partial mirror splits a beamlet into a
transmitted and a reflected child, so the cavity becomes a beam tree, and a [`Detector`](@ref) sums
the fields of all its leaves coherently.

## A Fabry-Pérot etalon

The etalon is a 10 µm thick glass plate whose faces carry lossless mirrors of 80 % reflectance. A flat
lens with infinite radii serves as the plate; its faces are coated with [`surface_model!`](@ref).

```@example fabry_perot
using GLMakie, BeamletOptics
const mm = 1e-3
const μm = 1e-6

n, d, R = 1.5, 10μm, 0.8
etalon = SphericalLens(Inf, Inf, d, 2mm, ConstantRefractiveIndex(n))
surface_model!(etalon, FrontFace(), BeamsplitterCoating(R))
surface_model!(etalon, BackFace(), BeamsplitterCoating(R))

detector = Detector(1mm)
translate3d!(detector, [0, 1mm, 0])
system = System([etalon, detector])
nothing # hide
```

The reflections never end. `depth_max` of [`solve_system!`](@ref) bounds the number of splits along a
branch of the beam tree: 40 splits are 20 round trips, after which `R⁴⁰ ≈ 10⁻⁴` of the field is left in
the cavity. We scan the wavelength of an [`AstigmaticGaussianBeamlet`](@ref) of 1 W and read the
transmitted power from the detector.

```@example fabry_perot
λs = range(615e-9, 635e-9; length = 201)
transmission = map(λs) do λ
    empty!(detector)
    beamlet = AstigmaticGaussianBeamlet([0, -1mm, 0], [0, 1, 0], λ, 100μm, 100μm; P0 = 1.0)
    solve_system!(system, beamlet; depth_max = 40)
    optical_power(detector)
end
nothing # hide
```

The transmission of an etalon follows the Airy function
``T = 1 / (1 + F \sin^2(2\pi n d / \lambda))`` with the coefficient of finesse ``F = 4R/(1-R)^2``
[BornWolf:1999](@cite):

```@example fabry_perot
airy(λ) = 1 / (1 + 4R / (1 - R)^2 * sin(2π * n * d / λ)^2)

fig = Figure(size = (650, 380))
ax = Axis(fig[1, 1]; xlabel = "wavelength [nm]", ylabel = "transmission", title = "Fabry-Pérot etalon")
lines!(ax, λs / 1e-9, airy.(λs); color = :black, label = "Airy function")
scatter!(ax, λs[1:4:end] / 1e-9, transmission[1:4:end]; color = :red, markersize = 7, label = "traced")
axislegend(ax; position = :rt)
fig
```

At 625 nm a round trip through the glass, ``2 n d = 30`` µm, is 48 wavelengths: all parts leave in phase
and the etalon transmits all the light, though each mirror reflects 80 %.

## A Fizeau wedge

Two plates with mirrors on their inner faces, 20 µm apart and tilted against each other by 0.1°, form a
wedge. Its thickness grows by ``\alpha`` per length, so the transmission repeats every
``\Delta x = \lambda / (2\alpha)``: the fringes of equal thickness of a Fizeau interferometer. Their
sharpness comes from the many reflections between the mirrors; with weak mirrors they would be the
sinusoidal fringes of two-beam interference.

```@example fabry_perot
α, λ = deg2rad(0.1), 632.8e-9

function mirror_plate(inner)
    plate = RectangularCompensatorPlate(5mm, 5mm, 0.5mm, ConstantRefractiveIndex(1.5))
    surface_model!(plate, FlatFace([0, inner, 0]), BeamsplitterCoating(R))
    surface_model!(plate, FlatFace([0, -inner, 0]), AntireflectionCoating())
    return plate
end

lower, upper = mirror_plate(1), mirror_plate(-1)
translate3d!(lower, [0, -0.5mm, 0])
zrotate3d!(lower, α / 2)
translate3d!(upper, [0, 20μm, 0])
zrotate3d!(upper, -α / 2)

screen = Detector(1.5mm)
translate3d!(screen, [0, 2mm, 0])

beamlet = AstigmaticGaussianBeamlet([0, -2mm, 0], [0, 1, 0], λ, 300μm, 300μm; P0 = 1.0)
solve_system!(System([lower, upper, screen]), beamlet; depth_max = 40)

x, z, I = intensity(screen, RegularGrid(num_points = 300))
fig = Figure(size = (560, 500))
ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x [mm]", ylabel = "z [mm]",
    title = "Fizeau fringes, Δx = λ/(2α) = $(round(λ / (2α) / mm; digits = 3)) mm")
heatmap!(ax, x / mm, z / mm, I; colormap = :inferno)
fig
```

The fringes run parallel to the edge of the wedge. Within the Gaussian beam of 0.3 mm the transmission
rises sharply wherever the local thickness is a multiple of half a wavelength.
