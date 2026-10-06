# Camera and scene helpers, implemented in ext/RenderCamera.jl

"""
    get_view(ls)

Returns the current camera view matrix of an `LScene` `ls`, if a suitable backend is loaded.

Handy at the REPL to freeze a view found interactively: rotate the scene by hand, call
`get_view(ax)`, and paste the printed matrix into the script as a literal passed to
[`set_view`](@ref).

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
get_view(::Any) = throw(MissingBackendError())

"""
    set_view(ls, view::AbstractMatrix)
    set_view(ls, eye, lookat, up)

Sets the camera of an `LScene` `ls`, if a suitable backend is loaded, either directly from
a view matrix (e.g. one obtained via [`get_view`](@ref)) or from an eye/lookat/up triple.
See also [`look_at!`](@ref) to aim the camera at a known point instead of specifying the
triple directly.

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
set_view(::Any, ::AbstractMatrix) = throw(MissingBackendError())
set_view(::Any, ::Any, ::Any, ::Any) = throw(MissingBackendError())

"""
    set_orthographic(ls)

Switches the camera of an `LScene` `ls` to an orthographic projection, if a suitable
backend is loaded. If not, a [`MissingBackendError`](@ref) will be thrown.
"""
set_orthographic(::Any) = throw(MissingBackendError())

"""
    hide_axis(ls, hide::Bool=true)

Hides the axis markers of an `LScene` `ls`, if a suitable backend is loaded. Can be
toggled via `hide`. If no suitable backend is loaded, a [`MissingBackendError`](@ref)
will be thrown.
"""
hide_axis(::Any, ::Bool) = throw(MissingBackendError())

"""
    arrow!(ax, pos, dir; scale=1, kwargs...)

Draws a single 3D arrow from `pos` pointing along `dir` into `ax`, if a suitable backend
is loaded, scaled to a fixed on-screen length (independent of `dir`'s own norm) so it
stays legible next to CAD geometry. If no suitable backend is loaded, a
[`MissingBackendError`](@ref) will be thrown.
"""
arrow!(::Any, ::AbstractVector, ::AbstractVector; kwargs...) = throw(MissingBackendError())

"""
    look_at!(ax, target, offset; up = [0, 0, 1])

Aims the camera of `ax` at `target` from `target + offset`, if a suitable backend is
loaded. A deterministic replacement for manually orbiting the scene to find a viewpoint,
handy for reproducible close-up figures. If no suitable backend is loaded, a
[`MissingBackendError`](@ref) will be thrown.
"""
look_at!(::Any, ::AbstractVector, ::AbstractVector; kwargs...) = throw(MissingBackendError())

"""
    render_lcs!(ax, pos, lcs; scale = 10, show_labels = false)
    render_lcs!(ax, object; scale = 10, show_labels = false)

Draws the local coordinate system of an object (or of an explicit `pos`/orientation pair)
into `ax` as a red/green/yellow arrow triad, if a suitable backend is loaded. Useful to
make the reference frame of an imported CAD mesh visible in the scene. If no suitable
backend is loaded, a [`MissingBackendError`](@ref) will be thrown.
"""
render_lcs!(::Any, ::AbstractArray = zeros(3), ::AbstractMatrix = Matrix{Float64}(I, 3, 3); kwargs...) =
    throw(MissingBackendError())
render_lcs!(::Any, ::AbstractObject; kwargs...) = throw(MissingBackendError())
