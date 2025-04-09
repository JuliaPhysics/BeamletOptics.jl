# Mirrors

A common optical element with a straight-forward optical interaction. This kind of component is in general defined as a [`BeamletOptics.AbstractReflectiveOptic`](@ref). For a basic [`Ray`](@ref) the interaction is simply defined by the [`BeamletOptics.reflection3d`](@ref) function. A more complex algorithm is required when when a [`PolarizedRay`](@ref) interacts with a reflecting surface. The polarization calculus that is performed is explained in the [Polarized Rays](@ref) section. Below, some of the concrete implemented mirror types are shown. In general, the [`Mirror`](@ref) is used as a concrete type to represent an arbitrary reflecting shape.

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

```@eval
file_dir = joinpath(@__DIR__, "..", "assets")

Base.include(@__MODULE__, joinpath(file_dir, "plano_mirror_showcase.jl"))

take_screenshot("plano_mirror_showcase.png", system, beam; size=(600, 400), view=mirror_camera, color=RGBf(0,1,0), flen=.4)
```

![Plano mirror showcase](plano_mirror_showcase.png)

## Concave Mirrors

The [`ConcaveSphericalMirror`](@ref) represents an ideal optical element with a spherical concave reflective surface, commonly used for non-dispersive focusing applications. Its geometry is modeled using a combination of a concave spherical surface and a plano substrate, represented internally by a [`BeamletOptics.UnionSDF`](@ref) (refer also to the [SDF-based spherical lenses](@ref) section).

```@eval
file_dir = joinpath(@__DIR__, "..", "assets")

Base.include(@__MODULE__, joinpath(file_dir, "spherical_mirror_showcase.jl"))

save("concave_mirror_showcase.png", fig, px_per_unit=4); nothing
```

![Concave mirror multipass showcase](concave_mirror_showcase.png)

The following constructor allows the spawning of concave spherical mirrors.

```@docs; canonical=false
ConcaveSphericalMirror(::Real, ::Real, ::Real)
```

