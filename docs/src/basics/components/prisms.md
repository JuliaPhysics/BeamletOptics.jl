# Prisms

Refractive prisms are represented by the [`Prism`](@ref) type, which wraps any closed shape together with a (wavelength dependent) refractive index. The package provides constructors for prisms with a convex polygonal cross-section, i.e. right prisms: a 2D polygon in the local x-y-plane which is extruded symmetrically along the local z-axis (the `height`). The local origin is the origin of the polygon coordinates, and like all components the prisms can be moved freely, see [Moving optical elements](@ref).

!!! info "Alignment"
    Prisms are *not aligned* with the global optical axis (+y) after construction. Use the kinematic functions to orient them, e.g. `zrotate3d!` to rotate about the extrusion axis.

## Arbitrary convex prisms

[`Prism`](@ref) accepts the corner points of the cross-section in [m] in either orientation (they are normalized to counter-clockwise). The polygon must be strictly convex, otherwise an `ArgumentError` is thrown. The underlying shape is the [`BeamletOptics.PolygonPrismSDF`](@ref), whose signed distance function is exact. All faces are flat, so rays enter and leave exactly on the faces.

```julia
# a 30-60-90 prism, made from NBK7
prism = Prism([(0.0, 0.0), (30e-3, 0.0), (0.0, 17.32e-3)], 20e-3, NBK7)
```

```@docs; canonical=false
Prism(::AbstractVector, ::Real, ::BeamletOptics.RefractiveIndex)
BeamletOptics.PolygonPrismSDF
```

## Common prisms

The [`EquilateralPrism`](@ref) is centered at the centroid of its triangle with the apex along +y and is the classic dispersing prism. A ray traversing it symmetrically at minimum deviation is deflected by ``\delta = 2\arcsin(n \sin(\alpha/2)) - \alpha`` with the apex angle ``\alpha = 60^\circ``. The [`DovePrism`](@ref) is a trapezoid with 45° end faces whose long axis is the y-axis. The [`RightAnglePrism`](@ref) has its legs along x and y.

```@docs; canonical=false
EquilateralPrism(::Real, ::Real, ::BeamletOptics.RefractiveIndex)
DovePrism(::Real, ::Real, ::Real, ::BeamletOptics.RefractiveIndex)
RightAnglePrism(::Real, ::Real, ::BeamletOptics.RefractiveIndex)
```
