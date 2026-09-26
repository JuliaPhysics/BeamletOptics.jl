module TestViewCube

using BeamletOptics
using Makie
using Test
using LinearAlgebra: normalize, norm, dot

const BMO = BeamletOptics

@testset "View cube" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    # All 26 regions: 6 faces, 12 edges and 8 corners
    regions = [s for s in Iterators.product(-1:1, -1:1, -1:1) if s != (0, 0, 0)]
    nnz(s) = count(!iszero, s)

    @testset "ray/cube intersection" begin
        # hit from outside, nearest face
        p = Ext._ray_cube_hit([0.2, -0.3, 5.0], [0.0, 0.0, -1.0])
        @test p ≈ [0.2, -0.3, 1.0]
        p = Ext._ray_cube_hit([-5.0, 0.5, 0.5], normalize([1.0, 0.0, 0.0]))
        @test p ≈ [-1.0, 0.5, 0.5]
        # oblique hit on an edge
        p = Ext._ray_cube_hit([3.0, 0.0, 3.0], normalize([-1.0, 0.0, -1.0]))
        @test p ≈ [1.0, 0.0, 1.0]
        # miss: parallel outside, beside, pointing away
        @test isnothing(Ext._ray_cube_hit([0.0, 2.0, 5.0], [0.0, 0.0, -1.0]))
        @test isnothing(Ext._ray_cube_hit([0.0, 0.0, 5.0], normalize([1.0, 0.0, -0.1])))
        @test isnothing(Ext._ray_cube_hit([0.0, 0.0, 5.0], [0.0, 0.0, 1.0]))
        # from inside: the exit point
        p = Ext._ray_cube_hit([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        @test p ≈ [0.0, 1.0, 0.0]
        p = Ext._ray_cube_hit([0.5, 0.0, 0.0], [-1.0, 0.0, 0.0])
        @test p ≈ [-1.0, 0.0, 0.0]
    end

    @testset "region classification" begin
        @test length(regions) == 26
        @test count(s -> nnz(s) == 1, regions) == 6
        @test count(s -> nnz(s) == 2, regions) == 12
        @test count(s -> nnz(s) == 3, regions) == 8
        # the center of each region on the surface of the cube is the region itself
        for s in regions
            @test Ext._cube_region(Float64.(collect(s))) == s
        end
        # within the margin of a face, but close to an edge
        @test Ext._cube_region([0.65, 0.0, 1.0]) == (0, 0, 1)
        @test Ext._cube_region([0.75, 0.0, 1.0]) == (1, 0, 1)
        @test Ext._cube_region([-0.8, -0.9, -1.0]) == (-1, -1, -1)
    end

    @testset "view per region" begin
        faces = Dict((0, 0, 1) => [0, 1, 0], (0, 0, -1) => [0, 1, 0], (0, -1, 0) => [0, 0, 1],
            (0, 1, 0) => [0, 0, 1], (1, 0, 0) => [0, 0, 1], (-1, 0, 0) => [0, 0, 1])
        for (s, up) in faces
            o, u = Ext._region_view(s)
            @test o ≈ collect(s)
            @test u ≈ up
        end
        for s in regions
            o, u = Ext._region_view(s)
            @test o ≈ normalize(Float64.(collect(s)))
            @test norm(u) ≈ 1
            @test abs(dot(u, o)) < 1e-12
            # +z is up, unless looking along z
            nnz(s) > 1 && @test u[3] > 0
        end
    end

    @testset "rotation between directions" begin
        a, b = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]
        k, θ = Ext._rotation_between(a, b, [0.0, 0.0, 1.0])
        @test k ≈ [0, 0, 1]
        @test θ ≈ π / 2
        @test Ext._rotate(a, k, θ / 2) ≈ normalize([1.0, 1.0, 0.0])
        @test Ext._rotate(a, k, θ) ≈ b atol = 1e-12
        # antiparallel: about the fallback axis, made perpendicular to `a`
        k, θ = Ext._rotation_between(a, -a, [1.0, 0.0, 1.0])
        @test k ≈ [0, 0, 1]
        @test θ ≈ π
        @test Ext._rotate(a, k, θ) ≈ -a atol = 1e-12
        # fallback parallel to `a`: any perpendicular axis
        k, _ = Ext._rotation_between(a, -a, a)
        @test abs(dot(k, a)) < 1e-12 && norm(k) ≈ 1
    end

    function _scene()
        fig = Figure(; size = (600, 500))
        ax = LScene(fig[1, 1])
        mesh!(ax, Rect3f(Point3f(-1), Vec3f(2)))
        return fig, ax
    end

    _cam(ax) = cameracontrols(ax.scene)
    _dir(ax) = normalize(Vector{Float64}(_cam(ax).eyeposition[]) .- Vector{Float64}(_cam(ax).lookat[]))
    _press!(ax) = (events(ax.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.press))
    _release!(ax) = (events(ax.scene).mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release))
    _tick!(ax, dt) = (events(ax.scene).tick[] = Makie.Tick(Makie.RegularRenderTick, 0, 0.0, dt))
    # Pixel position of the point `p` of the cube
    _px(cube, p) = Tuple(Float64.(Makie.shift_project(cube.scene, Point3f(p...))))

    @testset "construction and corners" begin
        fig, ax = _scene()
        cube = view_cube!(ax)
        @test cube isa Ext.ViewCube
        @test sprint(show, cube) == "ViewCube(top_right, 110 px)"
        vp, cvp = ax.scene.viewport[], cube.scene.viewport[]
        @test widths(cvp) == Vec(110, 110)
        @test maximum(cvp) == maximum(vp) .- 10
        @test cube.scene in ax.scene.children
        # the cube plots are not part of the main scene and ignore its clip planes
        @test all(p -> !(p in ax.scene.plots), cube.faces)
        @test all(p -> isempty(p.clip_planes[]), [cube.faces; cube.highlight])
        close(cube)

        cube = view_cube!(ax; size = 80, corner = :bottom_left)
        cvp = cube.scene.viewport[]
        @test widths(cvp) == Vec(80, 80)
        @test minimum(cvp) == minimum(ax.scene.viewport[]) .+ 10
        close(cube)
        @test_throws ArgumentError view_cube!(ax; corner = :center)
        @test_throws ArgumentError view_cube!(ax; size = 0)
        @test_throws ArgumentError view_cube!(ax; duration = -1)
    end

    @testset "camera sync" begin
        fig, ax = _scene()
        cube = view_cube!(ax)
        @test normalize(cube.eye) ≈ _dir(ax) atol = 1e-6
        for eye in ([5.0, 0.3, 1.0], [-2.0, 4.0, -3.0], [0.1, -6.0, 0.5])
            _cam(ax).eyeposition[] = Vec3f(eye...)
            @test normalize(cube.eye) ≈ _dir(ax) atol = 1e-6
            @test norm(cube.eye) ≈ 4
        end
        # the up vector follows as well
        set_view(ax, [0, -5, 0], [0, 0, 0], [1, 0, 0])
        @test cube.up ≈ [1, 0, 0] atol = 1e-6
        close(cube)
    end

    @testset "click on the top face" begin
        fig, ax = _scene()
        cube = view_cube!(ax; duration = 0.3)
        set_view(ax, [3, -4, 2], [0.5, 0.2, 0.1], [0, 0, 1])
        lookat = Vector{Float64}(_cam(ax).lookat[])
        dist = norm(Vector{Float64}(_cam(ax).eyeposition[]) .- lookat)
        events(ax.scene).mouseposition[] = _px(cube, (0, 0, 1))
        @test Ext._region_at_cursor(cube) == (0, 0, 1)
        # hovering highlights the region
        @test cube.hovered == (0, 0, 1)
        @test cube.highlight.visible[]
        _press!(ax)
        _release!(ax)
        # animated, driven by ticks
        @test !isnothing(cube.anim)
        _tick!(ax, 0.1)
        d = _dir(ax)
        @test !(d ≈ [0, 0, 1])
        for _ in 1:3
            _tick!(ax, 0.1)
        end
        @test isnothing(cube.anim)
        d = _dir(ax)
        @test abs(d[1]) < 1e-6 && abs(d[2]) < 1e-6 && d[3] > 0
        @test Vector{Float64}(_cam(ax).upvector[]) ≈ [0, 1, 0] atol = 1e-6
        @test Vector{Float64}(_cam(ax).lookat[]) ≈ lookat atol = 1e-6
        @test norm(Vector{Float64}(_cam(ax).eyeposition[]) .- lookat) ≈ dist rtol = 1e-6
        # the cube follows the main camera
        @test normalize(cube.eye) ≈ [0, 0, 1] atol = 1e-6
        # leaving the cube removes the highlight
        events(ax.scene).mouseposition[] = (5.0, 5.0)
        @test isnothing(cube.hovered)
        @test !cube.highlight.visible[]
        close(cube)
    end

    @testset "all regions switch the view" begin
        fig, ax = _scene()
        cube = view_cube!(ax; duration = 0)
        for s in regions
            # a view from the side of the region, such that it is visible on the cube
            o, u = Ext._region_view(s)
            set_view(ax, 5 .* normalize(o .+ 0.2 .* u), [0, 0, 0], [0, 0, 1])
            events(ax.scene).mouseposition[] = _px(cube, s)
            @test Ext._region_at_cursor(cube) == s
            _press!(ax)
            _release!(ax)
            @test _dir(ax) ≈ o atol = 1e-6
            @test Vector{Float64}(_cam(ax).upvector[]) ≈ u atol = 1e-6
        end
        # opposite view, animated via an intermediate direction
        close(cube)
        cube = view_cube!(ax; duration = 0.2)
        set_view(ax, [0, 0, -5], [0, 0, 0], [0, 1, 0])
        Ext._set_region_view!(cube, (0, 0, 1))
        _tick!(ax, 0.1)
        d = _dir(ax)
        @test abs(d[3]) < 1e-6
        _tick!(ax, 0.1)
        @test _dir(ax) ≈ [0, 0, 1] atol = 1e-6
        close(cube)
    end

    @testset "event consumption" begin
        fig, ax = _scene()
        cube = view_cube!(ax; duration = 0)
        probe = Ref(0)
        l = on(events(ax.scene).mousebutton, priority = 1) do event
            probe[] += 1
            return Consume(false)
        end
        # click on the cube: press and release are consumed
        events(ax.scene).mouseposition[] = _px(cube, (0, 0, 1))
        _press!(ax)
        @test cube.pressed
        _release!(ax)
        @test probe[] == 0
        @test !cube.pressed
        # click in the cube viewport, but beside the cube
        vp = cube.scene.viewport[]
        events(ax.scene).mouseposition[] = Tuple(Float64.(minimum(vp) .+ 1))
        @test isnothing(Ext._region_at_cursor(cube))
        _press!(ax)
        _release!(ax)
        @test probe[] == 2
        # click outside of the cube
        events(ax.scene).mouseposition[] = (50.0, 50.0)
        _press!(ax)
        _release!(ax)
        @test probe[] == 4
        off(l)
        close(cube)
    end

    @testset "close removes all listeners" begin
        fig, ax = _scene()
        ev = events(ax.scene)
        cam = _cam(ax)
        obs = (ev.mousebutton, ev.mouseposition, ev.tick, cam.eyeposition, cam.lookat,
            cam.upvector, ax.scene.viewport)
        n0 = map(o -> length(Makie.listeners(o)), obs)
        cube = view_cube!(ax)
        @test map(o -> length(Makie.listeners(o)), obs) != n0
        close(cube)
        @test map(o -> length(Makie.listeners(o)), obs) == n0
        @test isempty(cube.listeners)
        @test !(cube.scene in ax.scene.children)
        # closing twice is fine
        close(cube)
    end

    @testset "live view" begin
        # Beam along +y, mirror at 45° reflects it along +x onto the detector
        m = RoundPlanoMirror(25e-3, 5e-3)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        pd = Detector(5e-3)
        zrotate3d!(pd, -π / 2)
        translate3d!(pd, [0.1, 0.1, 0])
        beam = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3)
        # Picking requires a screen and is replaced
        pick_plot = Ref{Any}(nothing)
        gui = live_view(System([m, pd]), beam; throttle = false, trace_budget = Inf,
            pick = ax -> (pick_plot[], 0))
        cube = gui.view_cube
        @test cube isa Ext.ViewCube
        @test cube.scene in gui.ax.scene.children
        cube.duration = 0
        ax = gui.ax
        scene = ax.scene
        _click!() = (_press!(ax); _release!(ax))
        beside = Tuple(Float64.(minimum(scene.viewport[]) .+ 20))

        # a click on the cube neither selects ...
        pick_plot[] = gui.controls.h.handles[1].plots[1]
        events(scene).mouseposition[] = _px(cube, (0, 0, 1))
        _click!()
        @test isnothing(gui.controls.selected[])
        @test _dir(ax) ≈ [0, 0, 1] atol = 1e-6
        # ... nor deselects a component
        events(scene).mouseposition[] = beside
        _click!()
        @test gui.controls.selected[] === m
        pick_plot[] = nothing
        # the edge between the top and the front face, visible from the top
        events(scene).mouseposition[] = _px(cube, (0, -0.9, 1))
        _click!()
        @test gui.controls.selected[] === m
        @test _dir(ax) ≈ normalize([0, -1, 1]) atol = 1e-6
        # a click beside the cube deselects as usual
        events(scene).mouseposition[] = beside
        _click!()
        @test isnothing(gui.controls.selected[])

        # the clip planes of the live view do not affect the cube
        Ext._add_clip_plane!(gui, [0, 0.1, 0], [0, 1, 0]; select = false)
        @test !isempty(scene.theme.clip_planes[])
        @test all(p -> isempty(p.clip_planes[]), [cube.faces; cube.highlight])

        # close removes the listeners of the cube
        cube_listeners = copy(cube.listeners)
        _registered(l) = any(p -> p[2] === l.f, Makie.listeners(l.observable))
        @test all(_registered, cube_listeners)
        close(gui)
        @test isempty(cube.listeners)
        @test !any(_registered, cube_listeners)
        @test !(cube.scene in scene.children)

        gui = live_view(System([m, pd]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.5e-3);
            view_cube = false)
        @test isnothing(gui.view_cube)
        @test isempty(gui.ax.scene.children)
        close(gui)
    end
end

end
