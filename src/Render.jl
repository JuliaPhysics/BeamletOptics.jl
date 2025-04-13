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

render!(::Any, ::Any, kwargs...) = throw(MissingBackendError())

# """
#     render_object_normals!(axis, mesh::Mesh; l=0.01)

# Helper function to visualize the `mesh` normals of an object. The normals are displayed with a standard `l`ength of 0.01
# """
# function render_object_normals!(axis, mesh::Mesh; l = 0.01)
#     for fID in 1:size(mesh.faces)[1]
#         nml = normal3d(mesh, fID)
#         pos = @view mesh.vertices[mesh.faces[fID, 1], :]
#         vec = pos + l * nml
#         _render_object_normal!(axis, pos, vec; color=:blue)
#     end
#     return nothing
# end
# _render_object_normal!(::Any, ::AbstractVector, ::AbstractVector; color = :blue) = nothing

