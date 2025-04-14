abstract type RenderException <: Exception end

message(e::RenderException) = e.msg
showerror(io::IO, e::RenderException) = print(io, message(e))

mutable struct MissingBackendError <: RenderException
    msg::String
    function MissingBackendError()
        msg = "It appears no suitable Makie backend is loaded in this session."
        return new(msg)
    end
end

"""
    render!(axis, thing; kwargs...)

Test docs for `render!`, FIXME
"""
render!(::Any, ::Any, kwargs...) = throw(MissingBackendError())
