module BeamletOpticsMakieExt

using BeamletOptics
import BeamletOptics: render!, RenderException, _RenderTypes, get_view, set_view, hide_axis

const BMO = BeamletOptics

using Makie: Axis3, LScene, mesh!, surface!, lines!, RGBf, RGBAf, scatter!
using GeometryBasics: Point2, Point3
using AbstractTrees: PreOrderDFS
using MarchingCubes: MC, march

const _RenderEnv = Union{
    Axis3,
    LScene,
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

render!(::A, ::_RenderTypes; kwargs...) where A<:Any = throw(InvalidAxisError(A))

render!(::_RenderEnv, ::T; kwargs...) where T<:_RenderTypes = throw(RenderNotImplementedError(T))

# include order dependant!
include("RenderBeam.jl")
include("RenderGaussian.jl")
include("RenderSDF.jl")
include("RenderMesh.jl")
include("RenderObjects.jl")
include("RenderLenses.jl")
include("RenderCylinderLenses.jl")
include("RenderPresets.jl")

"""
    get_view(ls::LScene)

Returns the current `eyeposition`, `lookat` and `upvector` of the scene `ls`.
"""
function get_view(ls::LScene)
    cam = ls.scene.camera_controls
    eye = cam.eyeposition[]
    lookat = cam.lookat[]
    up = cam.upvector[]
    return eye, lookat, up
end

"""
    set_view(ls::LScene, eye, lookat, up)

Sets the current `eyeposition`, `lookat` and `upvector` of the scene `ls`.
"""
function set_view(ls::LScene, eye, lookat, up)
    cam = ls.scene.camera_controls
    cam.eyeposition[] = eye
    cam.lookat[] = lookat
    cam.upvector[] = up
end

"""
    hide_axis(ls::LScene, hide::Bool=true)

Hides the axis markers in the `LScene`. Can be toggled via `hide`.
"""
hide_axis(ls::LScene, hide::Bool=true) = (ls.show_axis[] = !hide)

end