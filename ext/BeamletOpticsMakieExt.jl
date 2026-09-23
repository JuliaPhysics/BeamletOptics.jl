module BeamletOpticsMakieExt

using BeamletOptics
import BeamletOptics: render!, RenderException, _RenderTypes, get_view, set_view, hide_axis,
                       set_orthographic, arrow!, render_lcs!, look_at!,
                       AbstractRenderHandle, live_render!, update_render!, remove_render!,
                       pick_object, kinematic_controls!, live_view

const BMO = BeamletOptics

using Makie: Axis3, LScene, mesh!, surface!, lines!, RGBf, RGBAf, scatter!, text!,
             update_cam!, cameracontrols, arrows3d!
using GeometryBasics: Point2, Point3, Point3f, Vec3f
using AbstractTrees: PreOrderDFS
using MarchingCubes: MC, march
using LinearAlgebra: dot, cross, normalize, norm

const _RenderEnv = Union{
    Axis3,
    LScene
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
        if !(t <: _RenderTypes)
            throw(ErrorException("Type $t not supported"))
        end
        msg = "Render function not implemented for type $t"
        return new(msg, t)
    end
end

render!(::A, ::_RenderTypes; kwargs...) where {A <: Any} = throw(InvalidAxisError(A))

function render!(::_RenderEnv, ::T; kwargs...) where {T <: _RenderTypes}
    throw(RenderNotImplementedError(T))
end

# include order dependant!
include("RenderBeam.jl")
include("RenderPolarization.jl")
include("RenderGaussian.jl")
include("RenderAstigmaticGaussian.jl")
include("RenderSDF.jl")
include("RenderMesh.jl")
include("RenderObjects.jl")
include("RenderLenses.jl")
include("RenderCylinderLenses.jl")
include("RenderMirrors.jl")
include("RenderPresets.jl")
include("RenderPolarizers.jl")
include("RenderCamera.jl")
# live rendering, must come after all static renderers
include("LiveObjects.jl")
include("LiveBeams.jl")
include("LiveInteraction.jl")
include("LiveView.jl")

end
