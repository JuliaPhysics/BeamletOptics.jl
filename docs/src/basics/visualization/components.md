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

## Bounding spheres

A shape can state a sphere that encloses it, see [`BeamletOptics.bounding_sphere_of`](@ref). The solver then skips the shape for every ray that misses the sphere. The spheres are drawn as magenta wireframes, either along with the components via the keyword `show_bounding_sphere` of `render!`, or on their own via [`render_bounding_sphere!`](@ref), which also takes a single shape or a sphere value. This is a debugging aid, e.g. to check that the sphere of a custom shape encloses the shape and is tight. An object that consists of several parts, e.g. a doublet lens, and an object group show one sphere per part and, in orange (keyword `main_color`), the main sphere around all parts, which the solver tests first. A shape or object without a sphere (`NoBoundingSphere`), e.g. a `MeshDummy`, shows none.

```julia
using GLMakie, BeamletOptics

lens = SphericalLens(50e-3, -50e-3, 10e-3, 25.4e-3)
mirror = RoundPlanoMirror(25.4e-3, 5e-3)
translate3d!(mirror, [0, 50e-3, 0])
system = System([lens, mirror])

fig = Figure()
ax = LScene(fig[1, 1])
render!(ax, system; show_bounding_sphere = true)

# or separately, e.g. in another color
render!(ax, system)
render_bounding_sphere!(ax, system; color = :orange, linewidth = 2)
```

With [`live_render!`](@ref) and `show_bounding_sphere = true`, the spheres follow the components when they are moved, see [Live rendering](@ref).

```@docs; canonical=false
render_bounding_sphere!
```

## Shapes (advanced)

The following methods render the shapes that components are built from. They are mainly of interest when writing custom components, see the [Signed Distance Functions (SDFs)](@ref) and [Meshes](@ref) pages of the API documentation. A [`BeamletOptics.UnionSDF`](@ref) or [`BeamletOptics.DifferenceSDF`](@ref) is rendered as one mesh, which is merged from the analytical meshes of its SDFs if all of them have one, and otherwise sampled via marching cubes. The spherical, aspherical and acylindrical lens surfaces have dedicated analytical renderers that are used automatically when a lens is rendered.

```@docs
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractSDF)
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.ConicSDF)
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.DifferenceSDF{T, <:BeamletOptics.ConicSDF}) where T
render!(::Union{GLMakie.Axis3, GLMakie.LScene}, ::BeamletOptics.AbstractMesh)
```
