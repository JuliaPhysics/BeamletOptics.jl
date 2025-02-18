# Mirrors

A common optical element with a straight-forward optical interaction. This kind of component is in general defined as a [`SCDI.AbstractReflectiveOptic`](@ref). For a basic [`SCDI.Ray`](@ref) the interaction is simply defined by the [`SCDI.reflection3d`](@ref) function. A more complex algorithm is required when when a [`SCDI.PolarizedRay`](@ref) interacts with a reflecting surface. The polarization calculus that is performed is explained in the [Polarized Rays](@ref) section. Below, some of the concrete implemented mirror types are shown. In general, the [`SCDI.Mirror`](@ref) is used as a concrete type to represent an arbitrary reflecting shape.

```@docs; canonical=false
SCDI.Mirror
```

The following constructors can be used to generate flat reflecting shapes. Additional types are explained below.

- [`SCDI.SquarePlanoMirror2D`](@ref)
- [`SCDI.SquarePlanoMirror`](@ref)
- [`SCDI.RectangularPlanoMirror`](@ref)
- [`SCDI.Retroreflector`](@ref)


## Plano Mirrors

A category of mirrors with a flat reflecting surface. A round version of this mirror can be easily generated using the [`SCDI.RoundPlanoMirror`](@ref) or [`SCDI.RightAnglePrismMirror`](@ref) types:

```@docs; canonical=false
SCDI.RoundPlanoMirror(::Real, ::Real)
```

Below, a trivial example of a beam path propagating through a system of Ø1"-mirrors mounted in [KM100CP/M](https://www.thorlabs.de/thorproduct.cfm?partnumber=KM100CP/M#ad-image-0) kinematic mounts is shown (e.g. [PF10-03-P01](https://www.thorlabs.com/thorproduct.cfm?partnumber=PF10-03-P01)). Note that the mounts are modeled as [`SCDI.NonInteractableObject`](@ref)s.

```@eval
file_dir = joinpath(@__DIR__, "..", "assets")

Base.include(@__MODULE__, joinpath(file_dir, "plano_mirror_showcase.jl"))

save("plano_mirror_showcase.png", fig, px_per_unit=4); nothing
```

![Plano mirror showcase](plano_mirror_showcase.png)

## Concave Mirrors

The [`SCDI.ConcaveSphericalMirror`](@ref) represents an ideal optical element with a spherical concave reflective surface, commonly used for non-dispersive focusing applications. Its geometry is modeled using a combination of a concave spherical surface and a plano substrate, represented internally by a [`SCDI.UnionSDF`](@ref) (refer also to the [SDF-based spherical lenses](@ref) section).

```@eval
file_dir = joinpath(@__DIR__, "..", "assets")

Base.include(@__MODULE__, joinpath(file_dir, "spherical_mirror_showcase.jl"))

save("concave_mirror_showcase.png", fig, px_per_unit=4); nothing
```

![Concave mirror multipass showcase](concave_mirror_showcase.png)

The following constructor allows the spawning of concave spherical mirrors.

```@docs; canonical=false
SCDI.ConcaveSphericalMirror(::Real, ::Real, ::Real)
```

