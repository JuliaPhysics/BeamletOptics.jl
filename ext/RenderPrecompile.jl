# Precompiles the static and the live rendering of objects, systems and beams, such that the first
# `render!` or `live_render!` does not pay for the compilation. Runs without a Makie backend.
@setup_workload begin
    # Each beam type needs its own detector, since the hits of a detector are concretely typed
    function _render_precompile_system()
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        lens = SphericalLens(0.05, -0.05, 5e-3, 25.4e-3)
        translate3d!(lens, [0, 0.05, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        return System([m, lens, pd])
    end

    @compile_workload begin
        fig = Figure()
        ax = LScene(fig[1, 1])
        sys = _render_precompile_system()
        render!(ax, sys)
        h = live_render!(ax, sys)
        for beam in (Beam([0.0, 0, 0], [0.0, 1, 0]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3),
                CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40))
            solve_system!(_render_precompile_system(), beam)
            render!(ax, beam)
            hb = live_render!(ax, beam)
            update_render!(hb)
            remove_render!(hb)
        end
        translate3d!(rendered(first(render_children(h))), [0, 1e-3, 0])
        update_render!(h)
        remove_render!(h)
    end
end
