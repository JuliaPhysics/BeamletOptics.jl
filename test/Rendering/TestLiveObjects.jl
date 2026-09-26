module TestLiveObjects

using BeamletOptics
using Makie
using GeometryBasics
using Test
using LinearAlgebra

const BMO = BeamletOptics

# Applies the model matrix of a plot to the point p
function _apply_model(model::Makie.Mat4, p)
    q = model * Makie.Point4d(p[1], p[2], p[3], 1.0)
    return Makie.Point3d(q[1], q[2], q[3]) ./ q[4]
end

# Checks that the model matrix maps the reference pose P0, R0 onto the pose P, R
function _check_reference_points(plot, P0, R0, P, R; atol = 1e-4)
    got0 = _apply_model(Makie.transformation(plot).model[], P0)
    @test isapprox(collect(got0), collect(P); atol)
    for i in 1:3
        e = zeros(3)
        e[i] = 1.0
        ref = P0 + R0 * e
        want = P + R * e
        got = _apply_model(Makie.transformation(plot).model[], ref)
        @test isapprox(collect(got), collect(want); atol)
    end
end

@testset "Live rendering: objects & systems" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    fig = Figure()
    ax = LScene(fig[1, 1])

    @testset "SingleShape (SphericalLens, surface-based)" begin
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        h = live_render!(ax, lens)
        @test h isa Ext.ObjectRenderHandle
        @test length(h.plots) > 0
        n = length(h.plots)

        translate3d!(lens, [0.01, -0.02, 0.03])
        zrotate3d!(lens, deg2rad(12))
        xrotate3d!(lens, deg2rad(-7))
        update_render!(h)

        @test length(h.plots) == n # no regeneration
        P, R = BMO.position(lens), BMO.orientation(lens)
        for plot in h.plots
            _check_reference_points(plot, h.P0, h.R0, P, R)
        end

        # repeated updates: plot count and identities stay constant, no leaks
        ids_before = objectid.(h.plots)
        for k in 1:20
            translate3d!(lens, [1e-4, 0, 0])
            zrotate3d!(lens, deg2rad(0.5))
            update_render!(h)
        end
        @test length(h.plots) == n
        @test objectid.(h.plots) == ids_before

        remove_render!(h)
        @test isempty(h.plots)
    end

    @testset "SingleShape (RoundPlanoMirror, SDF marching-cubes mesh)" begin
        mir = RoundPlanoMirror(0.025, 0.005)
        h = live_render!(ax, mir)
        @test length(h.plots) == 1 + Ext._default_edges(:reflective) # mesh and the feature edges of the look
        @test h.plots[1] isa Makie.Mesh

        translate3d!(mir, [0.0, 0.02, 0.0])
        yrotate3d!(mir, deg2rad(30))
        update_render!(h)

        P, R = BMO.position(mir), BMO.orientation(mir)
        _check_reference_points(h.plots[1], h.P0, h.R0, P, R)

        remove_render!(h)
        @test isempty(h.plots)
    end

    @testset "Returning to the reference pose resets the transformation" begin
        mir = RoundPlanoMirror(0.025, 0.005)
        h = live_render!(ax, mir)
        translate3d!(mir, [0.5, 0.0, 0.0])
        update_render!(h)
        @test !(Makie.transformation(h.plots[1]).model[] ≈ Makie.Mat4d(I))
        translate3d!(mir, [-0.5, 0.0, 0.0])
        @test BMO.position(mir) == h.P0
        update_render!(h)
        @test Makie.transformation(h.plots[1]).model[] ≈ Makie.Mat4d(I)
        remove_render!(h)
    end

    @testset "Mesh-based object (RightAnglePrismMirror): full vertex check" begin
        # Compare with the rigid transformation, since a new render! would re-mesh the object
        prism = RightAnglePrismMirror(0.02, 0.01)
        h = live_render!(ax, prism; edges = false)
        @test length(h.plots) == 1
        plot = h.plots[1]

        # raw (untransformed) vertex data baked in at the reference pose
        raw_verts = copy(GeometryBasics.coordinates(plot[1][]))
        P0, R0 = h.P0, h.R0

        translate3d!(prism, [0.03, -0.01, 0.0])
        xrotate3d!(prism, deg2rad(25))
        update_render!(h)

        P, R = BMO.position(prism), BMO.orientation(prism)
        Rd = R * R0'
        model = Makie.transformation(plot).model[]
        for v in raw_verts
            got = _apply_model(model, v)
            want = Rd * (collect(Float64.(v)) - P0) + P
            @test isapprox(collect(got), want; atol = 1e-3)
        end
    end

    @testset "MultiShape (CubeBeamsplitter): rigid update + non-rigid fallback" begin
        cbs = CubeBeamsplitter(0.02, λ -> 1.5)
        h = live_render!(ax, cbs)
        n = length(h.plots)
        @test count(p -> p isa Makie.Mesh, h.plots) == 2 # merged front and back, coating
        @test n == 2 + Ext._default_edges(:refractive) # and the feature edges of the look

        translate3d!(cbs, [0.01, 0.02, -0.01])
        zrotate3d!(cbs, deg2rad(35))
        update_render!(h)
        @test length(h.plots) == n

        P, R = BMO.position(cbs), BMO.orientation(cbs)
        for plot in h.plots
            _check_reference_points(plot, h.P0, h.R0, P, R)
        end

        # break rigidity: move one sub-shape on its own
        old_ids = objectid.(h.plots)
        translate3d!(cbs.front, [0.05, 0.0, 0.0])
        update_render!(h)
        @test length(h.plots) == n
        @test objectid.(h.plots) != old_ids # fallback re-rendered the plots
        # reference pose was reset to the (post-fallback) current pose
        @test h.P0 == BMO.position(cbs)
        @test h.R0 == BMO.orientation(cbs)

        # moving another part than the first does not change the pose of the object
        old_ids = objectid.(h.plots)
        translate3d!(cbs.back, [0.0, 0.0, 0.01])
        update_render!(h)
        @test objectid.(h.plots) != old_ids

        remove_render!(h)
        @test isempty(h.plots)
    end

    @testset "pick_object" begin
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        h = live_render!(ax, lens)
        other = RoundPlanoMirror(0.025, 0.005)
        translate3d!(other, [0, 0.1, 0])
        h_other = live_render!(ax, other)

        @test pick_object(h, h.plots[1]) === lens
        @test pick_object(h, h_other.plots[1]) === nothing
        @test pick_object(h, nothing) === nothing

        # Picking returns primitive child plots of recipes
        child = lines!(h.plots[1], Makie.Point3f[(0, 0, 0), (1, 1, 1)])
        @test child.parent === h.plots[1]
        @test pick_object(h, child) === lens

        remove_render!(h)
        remove_render!(h_other)
    end

    @testset "AbstractObjectGroup: single rigid handle" begin
        group = ObjectGroup([RoundPlanoMirror(0.02, 0.004), RoundPlanoMirror(0.02, 0.004)])
        translate3d!(group.objects[2], [0, 0.05, 0])
        h = live_render!(ax, group)
        @test h isa Ext.ObjectRenderHandle
        n = length(h.plots)
        @test count(p -> p isa Makie.Mesh, h.plots) == 1 # one merged mesh
        @test n == 1 + Ext._default_edges(:reflective) # and the feature edges of the look

        translate3d!(group, [0.02, 0, 0])
        zrotate3d!(group, deg2rad(10))
        update_render!(h)
        @test length(h.plots) == n

        remove_render!(h)
    end

    @testset "System: one handle per object, update & pick" begin
        sys = System([RoundPlanoMirror(0.025, 0.005), SphericalLens(0.05, -0.05, 0.01, 0.02)])
        translate3d!(sys.objects[2], [0.0, 0.15, 0.0])
        hs = live_render!(ax, sys)
        @test hs isa Ext.SystemRenderHandle
        @test length(hs.handles) == 2

        translate3d!(sys.objects[1], [0.01, 0.0, 0.0])
        zrotate3d!(sys.objects[2], deg2rad(5))
        update_render!(hs)

        for (obj, oh) in zip(sys.objects, hs.handles)
            @test oh.obj === obj
            P, R = BMO.position(obj), BMO.orientation(obj)
            for plot in oh.plots
                _check_reference_points(plot, oh.P0, oh.R0, P, R)
            end
        end

        @test pick_object(hs, hs.handles[1].plots[1]) === sys.objects[1]
        @test pick_object(hs, hs.handles[2].plots[1]) === sys.objects[2]

        remove_render!(hs)
        @test all(isempty(oh.plots) for oh in hs.handles)
    end

    @testset "System with nested groups renders per leaf" begin
        m1 = RoundPlanoMirror(0.02, 0.004)
        m2 = RoundPlanoMirror(0.02, 0.004)
        m3 = RoundPlanoMirror(0.02, 0.004)
        translate3d!(m2, [0, 0.05, 0])
        translate3d!(m3, [0, 0.1, 0])
        inner = ObjectGroup([m2, m3])
        outer = ObjectGroup([m1, inner])
        sys = System([outer])
        h = live_render!(ax, sys)

        # one handle per leaf, hierarchy in parent
        @test length(h.handles) == 3
        @test [oh.obj for oh in h.handles] == [m1, m2, m3]
        @test all(oh -> oh isa Ext.ObjectRenderHandle, h.handles)
        @test h.parent[m1] === outer
        @test h.parent[m2] === inner
        @test h.parent[m3] === inner
        @test h.parent[inner] === outer
        @test !haskey(h.parent, outer)
        @test Ext._top_level(h, m2) === outer
        @test Ext._top_level(h, outer) === outer

        # pick_object returns the top-level object, _pick_leaf the rendered object
        plot_m3 = h.handles[3].plots[1]
        @test pick_object(h, plot_m3) === outer
        @test Ext._pick_leaf(h, plot_m3) === m3
        @test pick_object(h, nothing) === nothing

        # moving the outer group applies the transform to all leaf plots without re-rendering
        ids_before = [objectid.(oh.plots) for oh in h.handles]
        translate3d!(outer, [0.02, -0.01, 0.005])
        zrotate3d!(outer, deg2rad(15))
        update_render!(h)
        @test [objectid.(oh.plots) for oh in h.handles] == ids_before
        for oh in h.handles
            P, R = BMO.position(oh.obj), BMO.orientation(oh.obj)
            for plot in oh.plots
                _check_reference_points(plot, oh.P0, oh.R0, P, R)
            end
        end

        # moving a sub-object moves only its plots
        model_m1 = Makie.transformation(h.handles[1].plots[1]).model[]
        translate3d!(m2, [0, 0, 0.01])
        update_render!(h)
        @test [objectid.(oh.plots) for oh in h.handles] == ids_before
        @test Makie.transformation(h.handles[1].plots[1]).model[] == model_m1
        oh2 = h.handles[2]
        _check_reference_points(oh2.plots[1], oh2.P0, oh2.R0, BMO.position(m2), BMO.orientation(m2))

        remove_render!(h)
        @test all(isempty(oh.plots) for oh in h.handles)
    end

    @testset "pick_object after a fallback rerender" begin
        cbs = CubeBeamsplitter(0.02, λ -> 1.5)
        h = live_render!(ax, System([cbs]))
        old_plots = copy(h.handles[1].plots)
        translate3d!(cbs.front, [0.0, 0.0, 0.01])
        update_render!(h)
        @test all(p -> !any(q -> q === p, old_plots), h.handles[1].plots)
        @test all(p -> pick_object(h, p) === cbs, h.handles[1].plots)
        remove_render!(h)
    end

    @testset "show" begin
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        h = live_render!(ax, lens)
        @test occursin("ObjectRenderHandle", sprint(show, h))
        sys = System([RoundPlanoMirror(0.025, 0.005)])
        hs = live_render!(ax, sys)
        @test occursin("SystemRenderHandle", sprint(show, hs))
    end
end

end # module
