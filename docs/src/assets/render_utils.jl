get_view(ls::LScene) = ls.scene.camera.view[]
set_view(ls::LScene, new::AbstractMatrix) = (ls.scene.camera.view[] = new)

function set_orthographic(ax::LScene)
    cam = cameracontrols(ax.scene)
    cam.fov[] = 1
    return nothing
end

hide_axis(ls::LScene) = (ls.show_axis[] = false)

function take_screenshot(
    fname::String,
    system::BMO.Nullable{BMO.AbstractSystem},
    beam::BMO.Nullable{BMO.AbstractBeam};
    size::Tuple{Int, Int} = (600, 300),
    view::BMO.Nullable{AbstractMatrix} = nothing,
    flen::Real = 10e-2,
    px_per_unit::Int = 8,
    color::Union{Symbol, RGBf, RGBAf} = :red,
    optional_rendering_fct::BMO.Nullable{Function} = nothing
    )
    fig = Figure(; size)
    display(fig)
    ax = LScene(fig[1,1])
    if !isnothing(system)
        render!(ax, system)
    end
    if !isnothing(beam)
        render!(ax, beam; flen, color)
    end
    if !isnothing(view)
        set_view(ax, view)
    end
    if !isnothing(optional_rendering_fct)
        optional_rendering_fct(fig, ax)
    end
    hide_axis(ax)
    save(fname, fig; px_per_unit, update = false)
end

arrow!(ax::LScene, pos, dir; color=:blue) = arrows!(ax, [Point3(pos)], [Point3(dir*5e-3)], arrowsize=Vec3f(1e-3, 1e-3, 1e-3); color)