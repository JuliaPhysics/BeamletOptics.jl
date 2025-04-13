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

# include order dependant!
include("RenderBeam.jl")
include("RenderGaussian.jl")
include("RenderSDF.jl")
include("RenderMesh.jl")
include("RenderObjects.jl")
include("RenderPresets.jl")

function render_surface!(axis::_RenderEnv, X, Y, Z; kwargs...)
    surface!(axis, X, Y, Z; kwargs...)
end

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