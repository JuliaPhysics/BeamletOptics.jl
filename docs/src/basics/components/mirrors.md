```@setup mirrors
mirror_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "mirror_renders")

Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "plano_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "spherical_mirror_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "parabolic_mirror_showcase.jl"), use_placeholder=false)
Main.DocUtils.conditional_include(joinpath(mirror_showcase_dir, "oap_mirror_showcase.jl"), use_placeholder=false)
```

# Mirrors

A common optical element with a straight-forward optical interaction. This kind of component is in general defined as a [`BeamletOptics.AbstractReflectiveOptic`](@ref). For a basic [`Ray`](@ref) the interaction is simply defined by the [`BeamletOptics.reflection3d`](@ref) function. A more complex algorithm is required when when a [`PolarizedRay`](@ref) interacts with a reflecting surface. The polarization calculus that is performed is explained in the [Polarized rays](@ref) section. Below, some of the concrete implemented mirror types are shown. In general, the [`Mirror`](@ref) is used as a concrete type to represent an arbitrary reflecting shape.

```@docs; canonical=false
Mirror
```

The following constructors can be used to generate flat reflecting shapes. Additional types are explained below.

- [`SquarePlanoMirror2D`](@ref)
- [`SquarePlanoMirror`](@ref)
- [`RectangularPlanoMirror`](@ref)
- [`Retroreflector`](@ref)


## Plano Mirrors

A category of mirrors with a flat reflecting surface. A round version of this mirror can be easily generated using the [`RoundPlanoMirror`](@ref) or [`RightAnglePrismMirror`](@ref) types:

```@docs; canonical=false
RoundPlanoMirror(::Real, ::Real)
```

Below, a trivial example of a beam path propagating through a system of Ø1"-mirrors mounted in [KM100CP/M](https://www.thorlabs.de/thorproduct.cfm?partnumber=KM100CP/M#ad-image-0) kinematic mounts is shown (e.g. [PF10-03-P01](https://www.thorlabs.com/thorproduct.cfm?partnumber=PF10-03-P01)). Note that the mounts are modeled as [`NonInteractableObject`](@ref)s.

![Plano mirror showcase](plano_mirror_showcase.png)

## Spherical Mirrors

The [`SphericalMirror`](@ref) represents an ideal optical element with a spherical concave reflective surface, commonly used for non-dispersive focusing applications. Its geometry is modeled using a combination of a concave spherical surface and a plano substrate, represented internally by a [`BeamletOptics.UnionSDF`](@ref) (refer also to the [SDF-based spherical lenses](@ref) section).

![Spherical mirror multipass showcase](spherical_mirror_showcase.png)

The following constructor allows the spawning of spherical mirrors.

```@docs; canonical=false
SphericalMirror(::Real, ::Real, ::Real)
```

## Parabolic Mirrors

The [`ParabolicMirror`](@ref) represents an on-axis parabolic mirror. In contrast to the [`SphericalMirror`](@ref), a paraboloid focuses a collimated beam that is parallel to its optical axis into a single point without spherical aberration. Its surface $y = -\frac{x^2 + z^2}{4f}$ is the special case $x_{\text{off}} = 0$ of the [`BeamletOptics.OffAxisParaboloidSDF`](@ref), such that the focus lies at $(0, -f, 0)$.

![Parabolic mirror showcase](parabolic_mirror_showcase.png)

The following constructor allows the spawning of on-axis parabolic mirrors.

```@docs; canonical=false
ParabolicMirror(::Real, ::Real)
```

## Off-Axis Parabolic Mirrors

The [`OffAxisParabolicMirror`](@ref) represents an off-axis parabolic (OAP) mirror used for achromatic focusing and beam deflection without introducing spherical aberration.

Its geometry is constructed from a parent paraboloid with focal length $f$ and off-axis distance $x_{\text{off}}$, parameterized by the Reflected Focal Length ($RFL$) and deflection angle $\theta_d$ (default 90°):

$$f = RFL \cdot \cos^2\left(\frac{\theta_d}{2}\right), \quad x_{\text{off}} = RFL \cdot \sin(\theta_d)$$

![Off-Axis Parabolic mirror showcase](oap_mirror_showcase.png)

The following constructor allows the spawning of off-axis parabolic mirrors.

```@docs; canonical=false
OffAxisParabolicMirror(::Real, ::Real)
```
