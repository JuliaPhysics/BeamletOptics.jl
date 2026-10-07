# Bounding spheres of shapes, implemented in ext/RenderBoundingSphere.jl

"""
    render_bounding_sphere!(ax, sphere; color = :magenta, linewidth = 1, kwargs...)
    render_bounding_sphere!(ax, shape; color = :magenta, linewidth = 1, kwargs...)
    render_bounding_sphere!(ax, object; color = :magenta, main_color = :orange, linewidth = 1, kwargs...)
    render_bounding_sphere!(ax, system; color = :magenta, main_color = :orange, linewidth = 1, kwargs...)

Draws the bounding sphere of an `AbstractShape` into `ax` (an `LScene` or `Axis3`), if a suitable
backend is loaded: the sphere of [`BeamletOptics.bounding_sphere_of`](@ref) in world coordinates, for
the current position and orientation of the shape, as a wireframe of three great circles (parallel
to the xy-, yz- and zx-plane of the world frame) in one `lines` plot. The solver skips a shape for
every ray that misses this sphere, hence the function shows what the solver tests, e.g. to check the
`bounding_sphere_of` method of an own shape type: the sphere must enclose the whole shape and should
be tight.

- sphere: a [`SingleBoundingSphere`](@ref) or [`MultiBoundingSphere`](@ref) is drawn as it is, a
  [`NoBoundingSphere`](@ref) gives no plot
- shape: one plot, or none if the shape has no bounding sphere
- object: the sphere of its shape, or, if it consists of several parts (e.g. a doublet lens, a cube
  beamsplitter or an `ObjectGroup`), the sphere of each part, since the solver tests each part on its
  own, and the main sphere around all parts in `main_color`. The solver tests a ray against the main
  sphere first. An object without a main sphere, because a part has none, gets none. An object
  without a sphere, e.g. a `MeshDummy`, gets no plot.
- system: the spheres of all its objects

# Keyword args

- `color = :magenta`: color of the lines
- `main_color = :orange`: color of the main spheres of objects with several parts and of groups,
  also those nested in a group. Not passed to the plot.
- `linewidth = 1`: line width in screen units

All other `kwargs` are passed on to the `lines` plot of `Makie`. Returns `nothing`. The spheres are
drawn once and do not follow an object that is moved afterwards. For spheres that follow, use the
`show_bounding_sphere` option of [`render!`](@ref) for objects with [`live_render!`](@ref).

```julia
render!(ax, system)
render_bounding_sphere!(ax, system)

# or in one call
render!(ax, system; show_bounding_sphere = true)
```

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
render_bounding_sphere!(::Any, ::Union{AbstractBoundingSphere, AbstractShape, AbstractObject, AbstractSystem}; kwargs...) =
    throw(MissingBackendError())
