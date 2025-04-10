module BeamletOpticsMakieExt

using BeamletOptics
import BeamletOptics: render!, RenderException

const BMO = BeamletOptics

using Makie: Axis3, LScene, mesh!, surface!, lines!, RGBAf, scatter!
using GeometryBasics: Point2, Point3
using AbstractTrees: PreOrderDFS
using MarchingCubes: MC, march

const _RenderEnv = Union{
    Axis3,
    LScene,
}

const _RenderTypes = Union{
    BMO.AbstractRay,
    BMO.AbstractBeam,
    BMO.AbstractShape,
    BMO.AbstractObject,
    BMO.AbstractObjectGroup,
    BMO.AbstractSystem,
}

struct InvalidAxisError <: RenderException
    msg::String
    ax::Type
    function InvalidAxisError(ax::Type)
        msg = "Invalid axis input of type $ax, must be LScene or Axis3"
        return new(msg, ax)
    end
end

struct RenderNotImplementedError <: RenderException
    msg::String
    t::Type
    function RenderNotImplementedError(t::Type)
        if !(t<:_RenderTypes)
            throw(ErrorException("Type $t not supported"))
        end
        msg = "Render function not implemented for type $t"
        return new(msg, t)
    end
end

render!(::A, ::Any; kwargs...) where A<:Any = throw(InvalidAxisError(A))

render!(::_RenderEnv, ::T; kwargs...) where T<:Any = throw(RenderNotImplementedError(T))

include("RenderSDF.jl")
include("RenderMesh.jl")
include("RenderObjects.jl")
include("RenderPresets.jl")

# function _render_ray!(axis::_RenderEnv,
#         ray::AbstractRay,
#         ray_end::AbstractVector;
#         color = :blue,
#         show_pos = false)
#     lines!(axis,
#         [ray.pos[1], ray_end[1]],
#         [ray.pos[2], ray_end[2]],
#         [ray.pos[3], ray_end[3]],
#         color=color,
#         linewidth=1.0,
#         transparency=true)
#     if show_pos
#         scatter!(axis, ray.pos; color)
#     end
#     return nothing
# end

# function _render_beam!(axis, beam::Beam{T, R}; color=:red, flen=1.0) where {T<:Real, R<:PolarizedRay{T}}
#     I0 = intensity(first(rays(beam)))
#     for child in PreOrderDFS(beam)
#         for ray in rays(child)
#             I = intensity(ray)
#             render_ray!(axis, ray, color = RGBAf(1,0,0, I/I0), flen = flen)
#         end
#     end
#     return nothing
# end

# function render_surface!(axis::_RenderEnv, X, Y, Z; kwargs...)
#     surface!(axis, X, Y, Z; kwargs...)
# end

# function _render_object_normal!(axis::_RenderEnv,
#         pos::AbstractVector,
#         vec::AbstractVector;
#         color = :blue)
#     lines!(axis,
#         [pos[1], vec[1]],
#         [pos[2], vec[2]],
#         [pos[3], vec[3]],
#         color=color)
# end

end
