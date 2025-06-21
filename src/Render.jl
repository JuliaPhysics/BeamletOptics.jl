abstract type RenderException <: Exception end

message(e::RenderException) = e.msg
showerror(io::IO, e::RenderException) = print(io, message(e))

"""
    MissingBackendError

Custom `Exception` type that indicates that the `Makie` extension of `BeamletOptics` has not been loaded correctly.
"""
mutable struct MissingBackendError <: RenderException
    msg::String
    function MissingBackendError()
        msg = "It appears no suitable Makie backend is loaded in this session."
        return new(msg)
    end
end

"""A collection of all types from BMO which might be renderable in principle."""
const _RenderTypes = Union{
    AbstractRay,
    AbstractBeam,
    AbstractShape,
    AbstractObject,
    AbstractObjectGroup,
    AbstractSystem,
}

"""
    render!(axis, thing; kwargs...)

The `render!` function allows for the visualization of optical system and beams under the condition that a suitable backend is loaded.
This means that either one of the following packages must be loaded in combination with `BeamletOptics` via `using`:

1. GLMakie 
    - preferred for 3D viewing
    - use `LScene` or `Axis3` environments
2. CairoMakie 
    - preferred for the generation of high-quality .pngs
    - only `Axis3` is supported

If no suitable backend is loaded, a [`MissingBackendError`](@ref) will be thrown.

# Implementations reqs.

All concrete implementations of `render!` must adhere to the following minimal interface:

`render!(axis, thing; kwargs...)`

- `axis`: an axis type of the union of `LScene` or `Axis3`
- `thing`: an abstract or concrete object or beam type
- `kwargs`: custom or `Makie` keyword arguments that are passed to the underlying backend

Refer to the `BeamletOptics` extension docs for `Makie` for more information.
"""
render!(::Any, ::_RenderTypes, kwargs...) = throw(MissingBackendError())