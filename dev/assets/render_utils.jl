get_view(ls::LScene) = ls.scene.camera.view[]
set_view(ls::LScene, new::AbstractMatrix) = (ls.scene.camera.view[] = new)

function set_orthographic(ax::LScene)
    cam = cameracontrols(ax.scene)
    cam.fov[] = 1
    return nothing
end

hide_axis(ls::LScene) = (ls.show_axis[] = false)

arrow!(
    ax::LScene,
    pos::AbstractVector,
    dir::AbstractVector;
    # Custom kwargs
    scale=1,
    # Makie kwargs
    color=:blue,
    tiplength=0.2,
    tipradius=0.1,
    ) = arrows3d!(
        ax,
        [Point3(pos)],
        [Point3(dir*5e-3*scale)];
        tiplength,
        tipradius,
        color
)

function render_sphere!(ax, pos, r; kwargs...)
    # catch small radii render bug by capping min. radius
    if r < 1e-5
        r = 1e-5
    end
    render!(ax, BMO.SphereSDF(Point3{Float64}(pos), Float64(r)); transparency=true, kwargs...)
end

lens_color() = RGBf(0.678, 0.847, 0.902)
lens_color(alpha) = RGBAf(0.678, 0.847, 0.902, alpha)