"""
    render!(ax, object; material = nothing, edges = nothing, kwargs...)

Renders the `object` into the specified `ax`is. Additional `kwargs` can be piped through to the backend.

# Look

Each object is rendered with the material preset of its component class, e.g. clear glass for
lenses and prisms or silver for mirrors. The parts of composite objects (e.g. the prisms and the
coating of a `CubeBeamsplitter`) are rendered with the class of each part. Explicit kwargs like
`color`, `alpha` or `transparency` override the preset, for all parts of the object.

- `material = nothing`: one of `:refractive`, `:reflective`, `:coating`, `:polarizer`,
  `:detector`, `:mechanics` and `:interface`, overrides the component class
- `edges = nothing`: draws the feature edges of the object, i.e. the sharp edges and the boundary of
  open surfaces, as thin dark lines. By default in the `:cad` look for all materials except
  `:mechanics`, in the `:modern` look as faint silhouettes of the glass (`:refractive`, `:coating`
  and `:interface`) only, see [`set_render_look`](@ref)

The cemented interfaces of doublet and triplet lenses are rendered as thin amber surfaces. See also
[`studio_lighting!`](@ref) for the lighting of the scene.

# Examples

It is recommended to use the following snippet in order to generate plots:

```julia
using GLMakie, BeamletOptics

fig = Figure()
ax = LScene(fig[1,1]) # or Axis3
render!(ax, my_BMO_obj; color=:white)
```

Additional keyword arguments can be passed. Refer to the `Makie` and `BeamletOptics` documentation for supported options
for each `object`.
"""
render!(ax::_RenderEnv, object::BMO.AbstractObject; kwargs...) = _render!(ax, object; kwargs...)

# Dispatch helper fct. for RenderPolarizers.jl, do not remove
_render!(ax::_RenderEnv, obj::BMO.AbstractObject; kwargs...) = render!(ax, BMO.shape_trait_of(obj), obj; kwargs...)

"""
    render!(ax, ::SingleShape, obj; material = nothing, edges = nothing, kwargs...)

Renders the shape of `obj` as one mesh if it has an analytic tessellation (see `_has_mesh`),
otherwise via the `render!` method of the shape. The mesh plot gets the attributes of the
`material` (see `_material`), the `kwargs` override them. The feature `edges` are drawn for
analytic meshes only, by default depending on the material and the look, see `_default_edges`.
Mechanics given as a mesh (e.g. a `MeshDummy`) keep the smooth shading of their mesh. The
`show_normals` kwargs of `render!(ax, ::AbstractMesh)` also select the `render!` method of the shape.
"""
function render!(ax::_RenderEnv, ::BMO.SingleShape, obj; material = nothing,
        edges::Union{Nothing, Bool} = nothing, kwargs...)
    s = BMO.shape(obj)
    kw = (; _material(obj, material)..., kwargs...)
    # Mechanics (e.g. a housing from an STL file) have no feature edges by default, since the many
    # edges of a detailed mesh cover the optics, and keep the smooth shading of their mesh
    class = isnothing(material) ? _material_class(obj) : material
    mechanics = class === :mechanics
    edges = something(edges, _default_edges(class))
    if mechanics && s isa BMO.AbstractMesh
        render!(ax, s; kw...)
    elseif _has_mesh(s) && !haskey(kwargs, :show_normals) && !haskey(kwargs, :show_normals_length)
        _render_mesh!(ax, s; edges, cemented = obj isa BMO.Lens, kw...)
    else
        render!(ax, s; kw...)
    end
    return nothing
end

"""
    render!(ax, ::MultiShape, obj; edges = nothing, kwargs...)

Renders all parts of `obj`. The analytic meshes of the parts are merged into one mesh plot per set
of plot attributes (e.g. the material), see `_plot_collected!`, plus one plot of the feature
`edges` of the object. By default, each part contributes its edges according to its material
class and the look, see `_default_edges`, i.e. mechanics parts have no edges. Nested `MultiShape`
parts are merged into the meshes of the outermost object. The `kwargs`, e.g. `material` or
`color`, and an explicit `edges` apply to all parts.
"""
function render!(ax::_RenderEnv, ::BMO.MultiShape, obj; edges::Union{Nothing, Bool} = nothing,
        kwargs...)
    kw = isnothing(edges) ? kwargs : (; kwargs..., edges)
    if !isnothing(_MESH_COLLECTOR[])
        for _obj in BMO.shape(obj)
            render!(ax, _obj; kw...)
        end
        return nothing
    end
    parts = Any[]
    with(_MESH_COLLECTOR => parts) do
        for _obj in BMO.shape(obj)
            render!(ax, _obj; kw...)
        end
    end
    _plot_collected!(ax, parts)
    return nothing
end

"""
    render!(ax::_RenderEnv, sys::AbstractSystem)

Render all objects contained in the `sys`tem.
"""
function render!(ax::_RenderEnv, sys::BMO.AbstractSystem; kwargs...)
    # Avoid use of objects(sys)
    for _obj in sys.objects
        render!(ax, _obj; kwargs...)
    end
    return nothing
end
