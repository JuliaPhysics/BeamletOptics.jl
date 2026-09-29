# Rendering components and systems

Optical components are rendered by drawing their shape. Example renderings of the individual components can be found in the [Mirrors](@ref), [Lenses](@ref), [Beamsplitters](@ref) and [Polarizers](@ref) sections.

## Components

The generic method below covers every [`BeamletOptics.AbstractObject`](@ref), including groups of components such as an [`ObjectGroup`](@ref), whose members are rendered one after another.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractObject)
```

### Materials

The color, opacity and shading of a component follow its material class, e.g. clear glass for
lenses or metallic silver for mirrors, in the active look (`:modern` or `:cad`). The materials,
feature edges and lighting are described on the [Look](@ref) page. Each default can be overridden
with the corresponding keyword argument, e.g. `render!(ax, lens; color=:orange)`,
`render!(ax, mirror; transparency=true)` or `render!(ax, mount; material=:mechanics)`.

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
