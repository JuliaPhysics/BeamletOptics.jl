get_view(ls::LScene) = ls.scene.camera.view[]
set_view(ls::LScene, new::AbstractMatrix) = (ls.scene.camera.view[] = new)

hide_axis(ls::LScene) = (ls.show_axis[] = false)

function take_screenshot(
    fname::String,
    system::BMO.Nullable{BMO.AbstractSystem},
    beam::BMO.Nullable{BMO.AbstractBeam};
    size::Tuple{Int, Int} = (600, 300),
    view::BMO.Nullable{AbstractMatrix} = nothing,
    flen::Real = 10e-2,
    px_per_unit::Int = 8,
    color::Union{Symbol, RGBf, RGBAf} = :red
    )
    fig = Figure(; size)
    display(fig)
    ax = LScene(fig[1,1])
    if !isnothing(system)
        render_system!(ax, system)
    end
    if !isnothing(beam)
        render_beam!(ax, beam; flen, color)
    end
    if !isnothing(view)
        set_view(ax, view)
    end
    hide_axis(ax)
    save(fname, fig; px_per_unit, update = false)
end