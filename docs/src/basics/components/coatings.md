# Coatings

A coating changes what happens at a boundary: how much light is reflected and transmitted, with which
phase and polarization. In BMO a coating is a surface model, see [Surface channels](@ref). Until faces
of components can carry coatings, a coating is a thin layer, a [`CoatingLayer`](@ref), that is placed
flush on the face of a substrate or stands free. The solver finds the layer and the substrate at the
same boundary by itself, see [Coincident hits](@ref).

```@docs; canonical=false
CoatingLayer
```

!!! info "Orientation"
    The normal of the shape of a `CoatingLayer` points away from the substrate. The layers of a
    [`ThinFilmCoating`](@ref) are listed from the substrate outward, and the reflection phase of the
    lossless coatings depends on the side, see [Reflection phase at a beamsplitter coating](@ref).

## Ideal and spectral coatings

The lossless coatings reflect a given fraction of the power. A [`BeamsplitterCoating`](@ref) is
carried by every thin beamsplitter; [`MirrorCoating`](@ref) and [`AntireflectionCoating`](@ref) are its
limits. A [`SpectralCoating`](@ref) reflects as a function of the wavelength and the angle of
incidence, e.g. a dichroic mirror or the curve of a datasheet, also separately for s and p.

```@example coatings
using BeamletOptics

# A longpass dichroic: reflects below 800 nm, transmits above
longpass = SpectralCoating((λ, θ) -> λ < 800e-9 ? 0.98 : 0.02)

# A datasheet curve, interpolated linearly between the points
datasheet = SpectralCoating([450e-9, 550e-9, 650e-9], [0.08, 0.02, 0.06])

dichroic = CoatingLayer(BeamletOptics.RectangularFlatMesh(25e-3, 25e-3), longpass)
nothing # hide
```

```@docs; canonical=false
BeamsplitterCoating
MirrorCoating
AntireflectionCoating
SpectralCoating
```

## Thin-film coatings

A [`ThinFilmCoating`](@ref) computes reflection and transmission from the layers of a known design,
coherently, for s and p, at any angle and from both sides; absorbing layers take the extinction
coefficient of their material, see [`AbsorbingRefractiveIndex`](@ref). The phase of the reflected and
transmitted light follows from the design as well.

```@example coatings
λ0 = 550e-9
TiO2, MgF2 = ConstantRefractiveIndex(2.35), ConstantRefractiveIndex(1.38)

# Single-layer anti-reflection coating: R ≈ 1.3 % on glass instead of 4.3 %
ar = ThinFilmCoating(QuarterWaveLayer(MgF2, λ0); splitting = BeamletOptics.Dominant())

# Quarter-wave mirror (HL)⁸H, from the substrate outward: R > 99.9 % at λ0
hr = ThinFilmCoating(QuarterWaveLayer(TiO2, λ0),
    Iterators.flatten((QuarterWaveLayer(MgF2, λ0), QuarterWaveLayer(TiO2, λ0)) for _ in 1:8)...)
nothing # hide
```

By default a thin-film coating traces both the reflected and the transmitted light as child beams. For a
coating whose weak part is of no interest, e.g. an anti-reflection coating, `splitting =
BeamletOptics.Dominant()` traces the stronger part only.

```@docs; canonical=false
ThinFilmCoating
ThinFilmLayer
QuarterWaveLayer
```

Metal layers and other absorbing materials take a complex refractive index: the real index and its
extinction coefficient κ.

```@example coatings
silver = AbsorbingRefractiveIndex(ConstantRefractiveIndex(0.135), 3.985)   # at 633 nm
mirror = ThinFilmCoating(ThinFilmLayer(silver, 150e-9))
nothing # hide
```

```@docs; canonical=false
AbsorbingRefractiveIndex
extinction_coefficient
```

## Ghost reflections

An uncoated boundary refracts without tracing its weak Fresnel reflection. [`FresnelGhosts`](@ref)
traces both, e.g. to find the ghost images of a lens system; the keywords `power_cutoff` and `depth_max`
of [`solve_system!`](@ref) bound the growing beam tree.

```@docs; canonical=false
FresnelGhosts
```
