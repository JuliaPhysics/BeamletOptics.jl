```@setup prisms
prism_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "prism_assets")

Main.DocUtils.conditional_include(joinpath(prism_showcase_dir, "equilateral_prism_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(prism_showcase_dir, "dove_prism_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(prism_showcase_dir, "right_angle_prism_showcase.jl"))
```

# Prisms

Refractive prisms are represented by the [`Prism`](@ref) type, which wraps any closed shape together with a (wavelength dependent) refractive index. The package provides constructors for prisms with a convex polygonal cross-section, i.e. right prisms: a 2D polygon in the local x-y-plane which is extruded symmetrically along the local z-axis (the `height`). The local origin is the origin of the polygon coordinates, and like all components the prisms can be moved freely, see [Moving optical elements](@ref).

!!! info "Alignment"
    Except for the [`DovePrism`](@ref), whose long axis lies along +y after construction, prisms are *not aligned* with the global optical axis (+y). Use the kinematic functions to orient them, e.g. `zrotate3d!` to rotate about the extrusion axis.

## Arbitrary convex prisms

[`Prism`](@ref) accepts the corner points of the cross-section in [m] in either orientation (they are normalized to counter-clockwise). The polygon must be strictly convex, otherwise an `ArgumentError` is thrown. The underlying shape is the [`BeamletOptics.PolygonPrismSDF`](@ref), whose signed distance function is exact. All faces are flat, so rays enter and leave exactly on the faces.

```julia
# N-BK7, Schott Sellmeier coefficients
NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)

# a 30-60-90 prism
prism = Prism([(0.0, 0.0), (30e-3, 0.0), (0.0, 17.32e-3)], 20e-3, NBK7)
```

```@docs; canonical=false
Prism(::AbstractVector, ::Real, ::BeamletOptics.AbstractRefractiveIndex)
BeamletOptics.PolygonPrismSDF
```

## Common prisms

The following constructors create the prisms found in most catalogs. All figures show the prisms from above, i.e. looking along -z onto the local x-y-plane of the cross-section.

### Equilateral prism

The [`EquilateralPrism`](@ref) is centered at the centroid of its triangle with the apex along +y and is the classic dispersing prism. A ray traversing it symmetrically at minimum deviation is deflected by ``\delta = 2\arcsin(n \sin(\alpha/2)) - \alpha`` with the apex angle ``\alpha = 60^\circ``. With a dispersion model such as a [`SellmeierEquation`](@ref), ``n`` depends on the wavelength and the prism splits white light into its colors. The figure shows rays from 450 nm to 650 nm through an N-SF11 prism at minimum deviation for 550 nm: blue light is refracted the most.

![Equilateral prism dispersing white light](equilateral_prism_showcase.png)

```@docs; canonical=false
EquilateralPrism(::Real, ::Real, ::BeamletOptics.AbstractRefractiveIndex)
```

### Dove prism

The [`DovePrism`](@ref) is a trapezoid with 45° end faces whose long axis is the y-axis. A beam along the axis is refracted at the entrance face, totally internally reflected at the base and leaves parallel to its original direction. The image is mirrored: in the figure, the order of the three parallel rays is reversed after the prism. Rotating the prism about its long axis by an angle rotates the image by twice that angle.

![Dove prism reversing the order of three parallel rays](dove_prism_showcase.png)

```@docs; canonical=false
DovePrism(::Real, ::Real, ::Real, ::BeamletOptics.AbstractRefractiveIndex)
```

### Right angle prism

The [`RightAnglePrism`](@ref) has its legs along x and y. A ray that enters a leg at normal incidence is totally internally reflected at the hypotenuse and leaves through the other leg, deflected by 90°, as shown in the figure. A ray that enters through the hypotenuse is reflected at both legs and returned antiparallel.

![Right angle prism deflecting three parallel rays by 90°](right_angle_prism_showcase.png)

```@docs; canonical=false
RightAnglePrism(::Real, ::Real, ::BeamletOptics.AbstractRefractiveIndex)
```
