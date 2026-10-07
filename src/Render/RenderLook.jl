# Looks and lighting, implemented in ext/RenderLook.jl

"""
    look_colors() -> Dict{Symbol, RGBf}

The colors of the material classes (`:refractive`, `:reflective`, `:coating`, `:polarizer`,
`:detector`, `:mechanics`, `:interface`) of the active look, see [`set_render_look`](@ref), e.g. to
recolor the rendered objects of a class for a dark background.

Needs the `Makie` extension, i.e. a loaded Makie backend.
"""
function look_colors end

"""
    studio_lighting!(ax::Union{LScene, Axis3}; preset = :studio)

Sets up a CAD-like lighting rig in the 3D view `ax`, if a suitable backend is loaded: an ambient
light, a key light from the upper right front, a fill light from the left and a rim light from
behind, all relative to the camera. Backends with a single directional light (e.g. CairoMakie)
get the ambient and the key light only. `preset = :none` leaves the lights unchanged.
The GUI of the package BeamletOpticsGUI applies the rig by default, scenes created via [`render!`](@ref) call it
explicitly.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
studio_lighting!(::Any; kwargs...) = throw(MissingBackendError())

"""
    set_render_look(look::Symbol)

Sets the look of all subsequently rendered objects, `:modern` (default) or `:cad`. The `:modern`
look renders clear glass with faint silhouettes, metallic mirrors and neutral mechanics without
edge lines, the `:cad` look saturated materials with feature edge lines. Explicit kwargs of [`render!`](@ref), e.g.
`color`, `material` or `edges`, override the look.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
set_render_look(::Any) = throw(MissingBackendError())
