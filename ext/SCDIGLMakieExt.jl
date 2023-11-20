module SCDIGLMakieExt

using SCDI: faces, vertices, AbstractMesh, AbstractRay
import SCDI: render_object!,
    _render_ray!, render_gaussian_beam_surface!, _render_object_normal!, render_sdf_mesh!
using GLMakie: Axis3, mesh!, surface!, lines!
using GeometryBasics: Point2, Point3

function render_object!(axis::Axis3, mesh::AbstractMesh)
    mesh!(axis, vertices(mesh), faces(mesh), transparency = true)
    return nothing
end

function _render_ray!(axis::Axis3,
        ray::AbstractRay,
        ray_end::AbstractVector;
        color = :blue)
    lines!(axis,
        [ray.pos[1], ray_end[1]],
        [ray.pos[2], ray_end[2]],
        [ray.pos[3], ray_end[3]],
        color=color)
    return nothing
end

function render_gaussian_beam_surface!(axis::Axis3, X, Y, Z; kwargs...)
    surface!(axis, X, Y, Z; kwargs...)
end

function _render_object_normal!(axis::Axis3,
        pos::AbstractVector,
        vec::AbstractVector;
        color = :blue)
    lines!(axis,
        [pos[1], vec[1]],
        [pos[2], vec[2]],
        [pos[3], vec[3]],
        color=color)
end

render_sdf_mesh!(axis::Axis3, vertices, faces; transparency = true) = mesh!(axis, vertices, faces, transparency=transparency)

end
