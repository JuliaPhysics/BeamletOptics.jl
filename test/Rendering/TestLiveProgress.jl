module TestLiveProgress

using BeamletOptics
using Makie
using LinearAlgebra: normalize, cross, norm
using Test

const BMO = BeamletOptics

@testset "Live progress window" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    fig = Figure()
    ax = LScene(fig[1, 1])
    scatter!(ax, [Point3f(0), Point3f(1)])
    n0 = length(ax.scene.plots)

    @testset "construction" begin
        o = Ext._ProgressOverlay(ax)
        # In a child scene, drawn after all plots of the 3D scene
        @test length(ax.scene.plots) == n0
        @test o.hud in ax.scene.children && length(o.hud.plots) == length(o.plots) == 4
        # GLMakie sorts the plots by the z translation before drawing, which also puts them in front
        # of the 3D scene; `overdraw` would put transparent plots over them
        @test all(p -> Makie.transformationmatrix(p)[][3, 4] > Ext._PROGRESS_Z, o.plots)
        @test all(p -> Makie.transformationmatrix(p)[][3, 4] < Ext._CARD_Z, o.plots)
        @test all(p -> !p.visible[], o.plots)
        @test all(p -> !p.overdraw[], o.plots)
        @test all(p -> !p.inspectable[], o.plots)
        @test all(p -> !p.transparency[], o.plots)
        @test all(p -> isempty(p.clip_planes[]), o.plots)
        # A NaN anchor does not enter the scene limits
        @test all(isnan, o.anchor[])
    end

    @testset "show and hide" begin
        o = Ext._ProgressOverlay(ax)
        n = length(ax.scene.plots)
        track = Ext._PROGRESS_TRACK[1]
        left = Ext._PROGRESS_BAR_X

        Ext._show_progress!(o, [0.1, 0.2, 0.3], 0.5, "Tracing beams 50 %")
        @test length(ax.scene.plots) == n
        @test all(p -> p.visible[], o.plots)
        @test o.anchor[] == Ext._screen_anchor(ax.scene, [0.1, 0.2, 0.3])
        @test o.label[] == "Tracing beams 50 %"
        @test only(vcat(o.plots[4].text[])) == "Tracing beams 50 %"
        @test o.fill_size[][1] ≈ 0.5 * track
        # The fill stays left-aligned with the track
        @test o.fill_offset[][1] - o.fill_size[][1] / 2 ≈ left
        @test o.plots[3].markersize[] == o.fill_size[]

        for (fraction, expected) in ((-0.1, 0.0), (1.3, 1.0), (0.5, 0.5))
            Ext._show_progress!(o, [0.1, 0.2, 0.3], fraction, "x")
            @test o.fill_size[][1] ≈ expected * track
            @test o.fill_offset[][1] - o.fill_size[][1] / 2 ≈ left
        end

        # Anchor update, integer points are accepted
        Ext._show_progress!(o, [1, 2, 3], 0.5, "Detector field 50 %")
        @test o.anchor[] == Ext._screen_anchor(ax.scene, Point3f(1, 2, 3))
        @test all(p -> o.anchor[] in vcat(p[1][]), o.plots)
        @test o.label[] == "Detector field 50 %"

        # Unchanged values do not notify the plots
        count = Ref(0)
        on(_ -> count[] += 1, o.anchor)
        on(_ -> count[] += 1, o.fill_size)
        on(_ -> count[] += 1, o.label)
        Ext._show_progress!(o, [1, 2, 3], 0.5, "Detector field 50 %")
        @test count[] == 0

        Ext._hide_progress!(o)
        @test all(p -> !p.visible[], o.plots)
        @test length(ax.scene.plots) == n
        Ext._hide_progress!(o)
        @test all(p -> !p.visible[], o.plots)
        Ext._show_progress!(o, [1, 2, 3], 0.25, "y")
        @test all(p -> p.visible[], o.plots)
        @test length(ax.scene.plots) == n
    end

    @testset "window stays in the view" begin
        scene = ax.scene
        w, h = Makie.widths(scene.viewport[])
        cam = cameracontrols(scene)
        eye, look = Vector(cam.eyeposition[]), Vector(cam.lookat[])
        v = normalize(look - eye)
        right = normalize(cross(v, Vector(cam.upvector[])))
        gap, margin, panel = Ext._PROGRESS_GAP, Ext._PROGRESS_MARGIN, Ext._PROGRESS_PANEL
        # The whole panel (from `gap` beyond the anchor on) lies inside the view with the margin
        inside(a) = all(a .+ gap .>= margin - 1e-3) &&
                    all(a .+ gap .+ panel .<= Vec2f(w, h) .- margin .+ 1e-3)
        # A point in the view: its projection
        a = Ext._screen_anchor(scene, look)
        q = Makie.project(scene, :data, :pixel, Point3f(look))
        @test a ≈ Point2f(q[1], q[2]) && inside(a)
        # Far to the right: at the right edge
        a = Ext._screen_anchor(scene, look + 1e3 * norm(look - eye) * right)
        @test inside(a) && a[1] ≈ w - margin - gap - panel[1]
        # Behind the camera, to the right: at the right edge (its projection is mirrored)
        a = Ext._screen_anchor(scene, eye - v + 0.5 * right)
        @test inside(a) && a[1] ≈ w - margin - gap - panel[1]
        # No position (NaN): at the bottom edge
        a = Ext._screen_anchor(scene, Point3f(NaN))
        @test inside(a) && a[2] ≈ margin - gap
    end
end

end
