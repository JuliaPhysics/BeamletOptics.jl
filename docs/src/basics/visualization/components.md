# Rendering components and systems

Optical components are rendered by drawing their shape. Example renderings of the individual components can be found in the [Mirrors](@ref), [Lenses](@ref), [Beamsplitters](@ref) and [Polarizers](@ref) sections.

## Components

The generic method below covers every [`BeamletOptics.AbstractObject`](@ref), including groups of components such as an [`ObjectGroup`](@ref), whose members are rendered one after another.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractObject)
```

### Default colours

Several component types are rendered with a preset colour and transparency:

| Type | `color` | `transparency` | `alpha` |
|---|---|---|---|
| `AbstractRefractiveOptic` | `:white` | `true` | – |
| `AbstractReflectiveOptic` | `:silver` | `false` | – |
| `Lens`, `DoubletLens`, `TripletLens` | light blue `RGBf(0.678, 0.847, 0.902)` | `true` | `0.5` |
| `ThinBeamsplitter` | `:magenta` | `true` | – |
| `NonInteractableObject` | `:grey` | `false` | – |
| `IntersectableObject` | `:grey` | `true` | – |

Each default can be overridden with the corresponding keyword argument, e.g. `render!(ax, lens; color=:orange)` or `render!(ax, mirror; transparency=true)`.

### Polarizing filters

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::PolarizationFilter)
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::LinearPolarizer)
```

## Systems

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractSystem)
```

## Shapes (advanced)

The following methods render the shapes that components are built from. They are mainly of interest when writing custom components, see the [Signed Distance Functions (SDFs)](@ref) and [Meshes](@ref) pages of the API documentation. A [`BeamletOptics.UnionSDF`](@ref) is rendered by drawing each of its SDFs. The spherical, aspherical and acylindrical lens surfaces have dedicated analytical renderers that are used automatically when a lens is rendered.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractSDF)
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.ConicSDF)
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.DifferenceSDF{T, <:BeamletOptics.ConicSDF}) where T
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractMesh)
```
