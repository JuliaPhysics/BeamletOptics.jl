# Precompiles the live rendering and the interactive controls, such that the first `live_view` does
# not pay for the compilation. Runs without a Makie backend, hence picking is replaced.
@setup_workload begin
    # Each beam type needs its own detector, since the hits of a detector are concretely typed
    function _live_precompile_system()
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
        sys = _live_precompile_system()
        h = live_render!(ax, sys)
        for beam in (Beam([0.0, 0, 0], [0.0, 1, 0]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3),
                CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40))
            hb = live_render!(ax, beam)
            solve_system!(_live_precompile_system(), beam)
            update_render!(hb)
            remove_render!(hb)
        end

        # Select the mirror, then drag it and use the keys in both modes
        ctrl = kinematic_controls!(ax, h; throttle = false, pick = _ -> (h.handles[1].plots[1], 0))
        scene = ax.scene
        events(scene).mouseposition[] = (10.0, 10.0)
        for action in (Mouse.press, Mouse.release, Mouse.press)
            events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, action)
        end
        events(scene).mouseposition[] = (20.0, 15.0)
        events(scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        for key in (Keyboard.up, Keyboard.left, Keyboard.page_up, Keyboard.m, Keyboard.up,
                Keyboard.backspace, Keyboard.h, Keyboard.v, Keyboard.v)
            events(scene).keyboardbutton[] = Makie.KeyEvent(key, Keyboard.press)
        end
        events(scene).unicode_input[] = '+'
        close(ctrl)
        update_render!(h)

        gui = live_view(System([RoundPlanoMirror(25e-3, 5e-3), Detector(5e-3)]),
            GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3))
        close(gui)
    end
end
