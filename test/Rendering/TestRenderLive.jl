module TestRenderLive

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

# A movable marker whose outline does not select it, see `pickable_plots`
struct _Marker
    pos::Vector{Float64}
end
Base.position(x::_Marker) = x.pos
BMO.orientation(::_Marker) = Matrix{Float64}(I, 3, 3)
BMO.kinematic_trait_of(::_Marker) = BMO.Movable(BMO.Oriented())
# clicks on the outline reach what lies behind it
BMO.pickable_plots(::_Marker, plots) = filter(p -> !(p isa Makie.Lines), plots)

# An object without a shape, whose rendering throws
struct _Unrenderable <: BMO.AbstractObject{Float64} end

# An own system handle on the protocol, e.g. of a GUI that combines the handles of several systems
struct _Combined <: BMO.AbstractSystemRenderHandle
    sys::System
    children::Vector{BMO.AbstractObjectRenderHandle}
end
BMO.rendered(h::_Combined) = h.sys
BMO.render_children(h::_Combined) = h.children
BMO.render_parent(::_Combined, _) = nothing
Base.delete!(h::_Combined, oh::BMO.AbstractObjectRenderHandle) = (filter!(c -> c !== oh, h.children); h)

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
        @test BMO.render_parent(h, m2) === inner
        @test BMO.render_parent(h, inner) === outer
        @test BMO.render_parent(h, outer) === nothing

        # pick_object of the system returns the top-level object, of the child the rendered object
        plot_m3 = h.handles[3].plots[1]
        @test pick_object(h, plot_m3) === outer
        @test pick_object(h.handles[3], plot_m3) === m3
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

    @testset "Render handle protocol" begin
        m = RoundPlanoMirror(0.025, 0.005)
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        g = ObjectGroup([lens])
        sys = System([m, g])
        h = live_render!(ax, sys)
        @test h isa BMO.AbstractSystemRenderHandle
        @test BMO.rendered(h) === sys
        children = BMO.render_children(h)
        @test all(c -> c isa BMO.AbstractObjectRenderHandle, children)
        @test [BMO.rendered(c) for c in children] == [m, lens]
        @test BMO.render_parent(h, lens) === g
        @test BMO.render_parent(h, m) === nothing
        @test length(BMO.render_plots(h)) == sum(c -> length(BMO.render_plots(c)), children)

        # a thing drawn by a function follows its pose
        src = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 2e-3, 1e-6; num_rings = 2, num_rays = 40)
        hm = live_render!(ax, src) do
            scatter!(ax, [Makie.Point3f(position(src))]; color = :orange)
            lines!(ax, [Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0.01, 0)])
        end
        @test hm isa BMO.AbstractObjectRenderHandle
        @test BMO.rendered(hm) === src
        @test length(BMO.render_plots(hm)) == 2
        translate3d!(src, [0.01, 0, 0])
        update_render!(hm)
        P, R = hm.P, hm.R
        @test P ≈ position(src)
        _check_reference_points(BMO.render_plots(hm)[1], hm.P0, hm.R0, P, R)

        # adding and removing a child handle after live_render!
        push!(h, hm)
        @test BMO.render_children(h)[end] === hm
        @test pick_object(h, BMO.render_plots(hm)[1]) === src
        delete!(h, hm)
        @test !any(c -> c === hm, BMO.render_children(h))
        @test pick_object(h, BMO.render_plots(hm)[1]) === nothing
        @test !isempty(BMO.render_plots(hm)) # delete! keeps the plots
        remove_render!(hm)
        remove_render!(h)
    end

    @testset "Objects added to and removed from a system handle" begin
        in_scene(p) = any(q -> q === p, ax.scene.plots)
        m = RoundPlanoMirror(0.025, 0.005)
        sys = System([m])
        h = live_render!(ax, sys)

        # a single object
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        push!(sys, lens)
        new = live_render!(h, lens)
        @test length(new) == 1
        @test BMO.rendered(only(new)) === lens
        @test BMO.render_children(h)[end] === only(new)
        @test isnothing(BMO.render_parent(h, lens))
        @test pick_object(h, BMO.render_plots(only(new))[1]) === lens
        @test_throws ArgumentError live_render!(h, lens)
        @test length(BMO.render_children(h)) == 2

        # a group is rendered per object, with its hierarchy
        m2 = RoundPlanoMirror(0.02, 0.004)
        m3 = RoundPlanoMirror(0.02, 0.004)
        translate3d!(m3, [0, 0.05, 0])
        inner = ObjectGroup([m3])
        group = ObjectGroup([m2, inner])
        push!(sys, group)
        new = live_render!(h, group)
        @test [BMO.rendered(oh) for oh in new] == [m2, m3]
        @test length(BMO.render_children(h)) == 4
        @test BMO.render_parent(h, m2) === group
        @test BMO.render_parent(h, m3) === inner
        @test BMO.render_parent(h, inner) === group
        @test isnothing(BMO.render_parent(h, group))
        @test pick_object(h, BMO.render_plots(new[2])[1]) === group
        @test_throws ArgumentError live_render!(h, m3)
        # an object of a group can not be removed on its own
        @test_throws "ObjectGroup" remove_render!(h, m2)
        @test_throws ArgumentError remove_render!(h, inner)
        @test length(BMO.render_children(h)) == 4
        @test BMO.render_parent(h, m2) === group

        # the plots of the added objects follow them
        translate3d!(group, [0.02, -0.01, 0.005])
        zrotate3d!(group, deg2rad(15))
        update_render!(h)
        for oh in new
            P, R = BMO.position(oh.obj), BMO.orientation(oh.obj)
            for plot in oh.plots
                _check_reference_points(plot, oh.P0, oh.R0, P, R)
            end
        end

        # removing the group deletes its plots, handles and hierarchy
        plots = [p for oh in new for p in BMO.render_plots(oh)]
        @test all(in_scene, plots)
        delete!(sys, group)
        @test isnothing(remove_render!(h, group))
        @test [BMO.rendered(oh) for oh in BMO.render_children(h)] == [m, lens]
        @test !any(in_scene, plots)
        @test isnothing(BMO.render_parent(h, m2))
        @test isnothing(BMO.render_parent(h, m3))
        @test isnothing(BMO.render_parent(h, inner))
        # nothing happens for an object that is not rendered
        remove_render!(h, group)
        @test length(BMO.render_children(h)) == 2

        # a group that can not be rendered completely leaves the handle and the axis unchanged
        m4 = RoundPlanoMirror(0.02, 0.004)
        broken = ObjectGroup([m4, _Unrenderable()])
        n_plots = length(ax.scene.plots)
        @test_throws Exception live_render!(h, broken)
        @test length(ax.scene.plots) == n_plots
        @test length(BMO.render_children(h)) == 2
        @test isnothing(BMO.render_parent(h, m4))
        @test length(live_render!(h, m4)) == 1
        remove_render!(h, m4)

        # an own system handle removes objects through the protocol, but can not render them
        c = _Combined(sys, copy(BMO.render_children(h)))
        @test_throws ArgumentError live_render!(c, m4)
        lens_plots = copy(BMO.render_plots(BMO.render_children(h)[2]))
        remove_render!(c, lens)
        @test [BMO.rendered(oh) for oh in BMO.render_children(c)] == [m]
        @test !any(in_scene, lens_plots)
        remove_render!(h)
    end

    @testset "pickable_plots" begin
        x = _Marker([0.0, 0, 0])
        h = live_render!(ax, x) do
            scatter!(ax, [Makie.Point3f(0, 0, 0)])
            lines!(ax, [Makie.Point3f(0, 0, 0), Makie.Point3f(1, 0, 0)])
        end
        s, l = BMO.render_plots(h)
        @test pick_object(h, s) === x
        @test pick_object(h, l) === nothing
        remove_render!(h)
    end

    @testset "own system handle on the protocol" begin
        s1, s2 = System([RoundPlanoMirror(0.025, 0.005)]), System([RoundPlanoMirror(0.025, 0.005)])
        h1, h2 = live_render!(ax, s1), live_render!(ax, s2)
        c = _Combined(s1, [BMO.render_children(h1); BMO.render_children(h2)])
        m2 = only(s2.objects)
        translate3d!(m2, [0.01, 0, 0])
        update_render!(c)
        oh = only(BMO.render_children(h2))
        _check_reference_points(BMO.render_plots(oh)[1], oh.P0, oh.R0, position(m2), BMO.orientation(m2))
        @test pick_object(c, BMO.render_plots(oh)[1]) === m2
        remove_render!(c)
        @test isempty(BMO.render_plots(c))
    end

    @testset "render! of a system, a component and a beam into an LScene" begin
        fig2 = Figure()
        ax2 = LScene(fig2[1, 1])
        m = RoundPlanoMirror(0.025, 0.005)
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 0.1, 0])
        sys = System([m])
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
        solve_system!(sys, beam)
        lens = SphericalLens(0.05, -0.05, 0.01, 0.02)
        n = length(ax2.scene.plots)
        render!(ax2, sys)
        render!(ax2, lens)
        render!(ax2, beam; color = :red, flen = 0.05)
        @test length(ax2.scene.plots) > n + 2
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
