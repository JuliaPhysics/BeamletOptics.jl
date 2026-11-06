# Bounding spheres of shapes, implemented in ext/RenderBoundingSphere.jl

"""
    render!(ax, sphere::AbstractBoundingSphere; color = :magenta, linewidth = 1, kwargs...)

Draws the bounding `sphere` into `ax` (an `LScene` or `Axis3`), if a suitable backend is loaded, as
a wireframe of three great circles (parallel to the xy-, yz- and zx-plane of the world frame) in one
`lines` plot. A [`SingleBoundingSphere`](@ref) or [`MultiBoundingSphere`](@ref) is drawn as it is, a
[`NoBoundingSphere`](@ref) gives no plot.

The solver skips a shape for every ray that misses its sphere, see
[`BeamletOptics.bounding_sphere_of`](@ref), hence the plot shows what the solver tests, e.g. to check
the `bounding_sphere_of` method of an own shape type: the sphere must enclose the whole shape and
should be tight. The sphere is computed for the current pose of the shape and does not follow a shape
that is moved afterwards.

# Keyword args

- `color = :magenta`: color of the lines
- `linewidth = 1`: line width in screen units

All other `kwargs` are passed on to the `lines` plot of `Makie`. Returns `nothing`.

```julia
render!(ax, BeamletOptics.bounding_sphere_of(lens); color = :orange)

# the spheres of all objects, one per part and the main sphere of each object of several parts
render!(ax, system; show_bounding_sphere = true)
```

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
render!(::Any, ::AbstractBoundingSphere; kwargs...) = throw(MissingBackendError())
