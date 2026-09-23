# for docs refer to Render.jl

get_view(ls::LScene) = ls.scene.camera.view[]

set_view(ls::LScene, view::AbstractMatrix) = (ls.scene.camera.view[] = view)

function set_view(ls::LScene, eye, lookat, up)
    cam = ls.scene.camera_controls
    cam.eyeposition[] = Vec3f(eye...)
    cam.lookat[] = Vec3f(lookat...)
    cam.upvector[] = Vec3f(up...)
    update_cam!(ls.scene, cam)
    return nothing
end

function set_orthographic(ls::LScene)
    cam = cameracontrols(ls.scene)
    cam.fov[] = 1
    return nothing
end

hide_axis(ls::LScene, hide::Bool = true) = (ls.show_axis[] = !hide)

look_at!(ls::LScene, target::AbstractVector, offset::AbstractVector; up = [0, 0, 1]) =
    set_view(ls, target .+ offset, target, up)

arrow!(
    ax::LScene,
    pos::AbstractVector,
    dir::AbstractVector;
    # kwargs
    scale = 1,
    # Makie kwargs
    color = :blue,
    tiplength = 0.2,
    tipradius = 0.1,
    kwargs...
) = arrows3d!(
    ax,
    [Point3(pos)],
    [Point3(dir * 5e-3 * scale)];
    tiplength,
    tipradius,
    color,
    kwargs...
)

function render_lcs!(
        ax::LScene,
        LCS_pos::AbstractArray = Point3(0),
        LCS::AbstractMatrix = [1 0 0; 0 1 0; 0 0 1];
        scale::Real = 10,
        show_labels::Bool = false
)
    arrow!(ax, LCS_pos, LCS[:,1]; scale, color=:red)
    arrow!(ax, LCS_pos, LCS[:,2]; scale, color=:green)
    arrow!(ax, LCS_pos, LCS[:,3]; scale, color=:yellow)
    if show_labels
        text!(ax, LCS_pos .+ scale * Point3(5e-3, 0, 1e-3), text="x")
        text!(ax, LCS_pos .+ scale * Point3(0, 5e-3, 1e-3), text="y")
        text!(ax, LCS_pos .+ scale * Point3(0, 0, 6.5e-3), text="z")
    end
    return nothing
end

function render_lcs!(ax::LScene, object::BMO.AbstractObject; scale::Real = 10, show_labels::Bool = false)
    render_lcs!(ax, BMO.position(object), BMO.orientation(object); scale, show_labels)
    return nothing
end
