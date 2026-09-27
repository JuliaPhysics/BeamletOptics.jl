# Lenses & surfaces

**Orientation:** optical axis +y, front vertex at y = 0, back vertex at y = `thickness(lens)`.
**Radius sign (ISO 10110):** R > 0 if the center of curvature lies at +y, R < 0 if at −y, `Inf` = plano.

## Spherical catalog lenses

```julia
SphericalLens(r1, r2, l, d = BeamletOptics.inch, n = λ -> 1.5)   # l: center thickness, d: diameter
ThinLens(R1, R2, d, n)                                          # thin-lens Lens; SphericalLens uses it for l = 0
SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2)             # cemented doublet
SphericalTripletLens(r1, r2, r3, r4, l1, l2, l3, d, n1, n2, n3) # cemented triplet
thickness(lens)                                                 # total center thickness
```

Examples (Thorlabs, N-BK7, `inch = BeamletOptics.inch`):

```julia
LB1811 = SphericalLens(34.9mm, -34.9mm, 6.8mm, inch, NBK7)   # biconvex
LA1805 = SphericalLens(Inf, -15.5mm, 8.6mm, inch, NBK7)      # plano-convex, flat side first
LD1464 = SphericalLens(-52mm, 52mm, 3mm, inch, NBK7)         # biconcave
LE1234 = SphericalLens(-82.2mm, -32.1mm, 3.6mm, inch, NBK7)  # meniscus
```

If the center thickness is too small for the radii and diameter, the constructor throws
"cylinder section length of ≤ 0, use ThinLens instead". Check the catalog values, increase `l`,
reduce `d`, or use `ThinLens`.

## General lenses from surfaces

```julia
Lens(front, back, center_thickness, n)   # both rotationally symmetric, or both cylindrical
Lens(front, center_thickness, n)         # plano back side
```

Rotationally symmetric surfaces:

- `SphericalSurface(radius, diameter[, mechanical_diameter])`
- `CircularFlatSurface(diameter)`
- `EvenAsphericalSurface(radius, diameter, conic_constant, coefficients[, mechanical_diameter])`

Cylindrical surfaces (curved in one direction, extruded over `height`):

- `CylindricalSurface(radius, diameter, height[, mechanical_diameter])`
- `AcylindricalSurface(radius, diameter, height, conic_constant, coefficients[, mechanical_diameter])`
- `RectangularFlatSurface(size)` for a plano side; both sides must have the same height.

**Aspheric coefficients** are the even terms `Σ α[i] r^(2i)` in SI units: `α[1]` is the r² term,
`α[2]` the r⁴ term, and so on. Convert mm-based catalog values like this:

```julia
A = [0, A4 * 1e3^3, A6 * 1e3^5, A8 * 1e3^7, A10 * 1e3^9]
AL75150 = Lens(EvenAsphericalSurface(76.68mm, 75mm, -0.675, A), 15mm, λ -> 1.5006520430)
```

Aspheres are experimental. General lenses support only spherical meniscus shapes.

## Multi-element lenses from parts

```julia
DoubletLens(front::Lens, back::Lens)
TripletLens(front::Lens, middle::Lens, back::Lens)
```

Pre-translate the elements yourself so that the cemented surfaces touch. Gaps give wrong results, and
TIR at cemented interfaces is not modeled. To model zoom or multi-group objectives, nest `ObjectGroup`s
of lenses (see `docs/src/assets/examples/lens_groups.jl` in the package repo).
