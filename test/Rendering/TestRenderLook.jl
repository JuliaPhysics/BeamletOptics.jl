module TestRenderLook

using BeamletOptics
using Makie
using GeometryBasics
using Test
using LinearAlgebra: norm, cross, I

const BMO = BeamletOptics

const Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

"""Returns the new plots of `render!(ax, x; kwargs...)`."""
function rendered_plots(x; kwargs...)
    ax = LScene(Figure()[1, 1])
    n0 = length(ax.scene.plots)
    render!(ax, x; kwargs...)
    return ax.scene.plots[(n0 + 1):end]
end

mesh_plots(plots) = filter(p -> p isa Makie.Mesh, plots)
edge_plots(plots) = filter(p -> p isa Makie.Lines, plots)

"""Returns the line segments of the polylines of the edge plot `p`, see `_polylines`."""
function segments(p::Makie.Lines)
    pts = p[1][]
    return [(pts[k], pts[k + 1]) for k in 1:(length(pts) - 1) if !any(isnan, pts[k]) && !any(isnan, pts[k + 1])]
end

points(p::Makie.Mesh) = [Vector{Float64}(q) for q in GeometryBasics.coordinates(p[1][])]
faces(p::Makie.Mesh) = [convert.(Int, Tuple(f)) for f in GeometryBasics.faces(p[1][])]
mesh_size(pts) = maximum(maximum(p[i] for p in pts) - minimum(p[i] for p in pts) for i in 1:3)

"""Returns true if every edge is shared by exactly two faces, after welding the vertices within `tol`."""
function is_closed(pts, fcs, tol)
    ids = Ext._weld(pts, tol)
    edges = Dict{Tuple{Int, Int}, Int}()
    for f in fcs
        w = map(i -> ids[i], f)
        allunique(w) || return false
        for k in 1:3
            key = minmax(w[k], w[mod1(k + 1, 3)])
            edges[key] = get(edges, key, 0) + 1
        end
    end
    return all(==(2), values(edges))
end

# Applies the model matrix of a plot to the point p
function apply_model(plot, p)
    q = Makie.transformation(plot).model[] * Makie.Point4d(p[1], p[2], p[3], 1.0)
    return [q[1], q[2], q[3]] ./ q[4]
end

"""Checks the look attributes of the mesh plot `p` against the material `mat`."""
function check_material(p, mat; color = mat.color)
    @test p.color[] == Makie.to_color(color)
    @test p.alpha[] ≈ mat.alpha
    @test p.transparency[] == mat.transparency
    @test p.diffuse[] ≈ Vec3f(mat.diffuse)
    @test p.specular[] ≈ Vec3f(mat.specular)
    @test p.shininess[] ≈ mat.shininess
end

@testset "Render look" begin
    @test !isnothing(Ext)
    mats = Ext._MATERIALS

    @testset "materials" begin
        # values of the table of the plan
        @test mats[:refractive].color == RGBf(0.72, 0.85, 0.92)
        @test mats[:refractive].alpha ≈ 0.35
        @test mats[:refractive].shininess ≈ 96
        @test Set(keys(mats)) == Set([:refractive, :reflective, :coating, :polarizer, :detector,
            :mechanics, :interface])

        lens = SphericalLens(34.9e-3, -34.9e-3, 6.8e-3, 25.4e-3)
        p = first(rendered_plots(lens))
        @test p isa Makie.Mesh
        check_material(p, mats[:refractive])
        # explicit kwargs override the material only partially
        p = first(rendered_plots(lens; color = :red))
        check_material(p, mats[:refractive]; color = :red)
        p = first(rendered_plots(lens; material = :mechanics))
        check_material(p, mats[:mechanics])
        @test_throws ArgumentError rendered_plots(lens; material = :gold)

        # component classes
        for (obj, class) in (RoundPlanoMirror(25e-3, 5e-3) => :reflective,
                Retroreflector(25e-3) => :reflective,
                Prism(BMO.BoxSDF(10e-3, 10e-3, 10e-3), λ -> 1.5) => :refractive,
                ThinBeamsplitter(10e-3, 10e-3) => :coating,
                RoundPolarizationFilter(25e-3) => :polarizer,
                Detector(10e-3) => :detector,
                NonInteractableObject(BMO.CylinderSDF(5e-3, 2e-3)) => :mechanics,
                IntersectableObject(BMO.CylinderSDF(5e-3, 2e-3)) => :mechanics)
            @test Ext._material_class(obj) == class
            check_material(first(rendered_plots(obj)), mats[class])
        end

        # LinearPolarizer: substrates refractive, film with the polarizer material, edges switchable
        lipo = RoundLinearPolarizer(25e-3, 2e-3, 2e-3, λ -> 1.5)
        meshes = mesh_plots(rendered_plots(lipo))
        @test count(p -> p.color[] == Makie.to_color(mats[:polarizer].color), meshes) == 1
        @test count(p -> p.color[] == Makie.to_color(mats[:refractive].color), meshes) == 2
        n_lines = length(edge_plots(rendered_plots(lipo)))
        n_lines_noedges = length(edge_plots(rendered_plots(lipo; edges = false)))
        @test n_lines - n_lines_noedges == 3 # one edge plot per part, the axis lines stay

        # MultiShape: materials per part, explicit kwargs for all parts
        cbs = CubeBeamsplitter(20e-3, λ -> 1.5)
        prisms, coating = mesh_plots(rendered_plots(cbs))
        check_material(prisms, mats[:refractive])
        check_material(coating, mats[:coating])
        prisms, coating = mesh_plots(rendered_plots(cbs; color = :blue))
        @test prisms.color[] == coating.color[] == Makie.to_color(:blue)
        @test coating.alpha[] ≈ mats[:coating].alpha
        prisms, coating = mesh_plots(rendered_plots(cbs; material = :reflective))
        check_material(prisms, mats[:reflective])
        check_material(coating, mats[:reflective])
    end

    @testset "feature edges" begin
        box = Prism(BMO.BoxSDF(20e-3, 10e-3, 30e-3), λ -> 1.5)
        zrotate3d!(box, deg2rad(35))
        xrotate3d!(box, deg2rad(-20))
        plots = rendered_plots(box)
        @test length(plots) == 2
        @test plots[1] isa Makie.Mesh
        edges = only(edge_plots(plots))
        @test length(segments(edges)) == 12
        @test edges.color[] == Ext._EDGE_COLOR
        @test edges.linewidth[] == 1
        @test length(Ext._feature_edges(Ext._tessellate(BMO.shape(box)))) == 2 * 12
        # the corners of the box end the polylines
        @test count(p -> any(isnan, p), edges[1][]) == 11

        # two circles, no edges along the smooth side wall
        s = BMO.CylinderSDF(5e-3, 2e-3)
        cyl = NonInteractableObject(s)
        edges = only(edge_plots(rendered_plots(cyl)))
        segs = segments(edges)
        @test length(segs) == 2 * Ext._N_THETA
        # both ends of each segment lie on the same end face (local y = ±height)
        y(p) = BMO._world_to_sdf(s, Vector{Float64}(p))[2]
        @test all(abs(abs(y(p)) - s.height) < 1e-8 for seg in segs for p in seg)
        @test all(sign(y(p)) == sign(y(q)) for (p, q) in segs)
        # two closed loops
        pts = edges[1][]
        @test count(p -> any(isnan, p), pts) == 1
        k = findfirst(p -> any(isnan, p), pts)
        @test pts[1] == pts[k - 1] && pts[k + 1] == pts[end]

        # zero thickness: boundary edges
        @test length(segments(only(edge_plots(rendered_plots(Detector(10e-3)))))) == 4

        # coincident flat faces (the bases of the caps of a thin lens) have no edges, only the rim
        tl = ThinLens(34.9e-3, -34.9e-3, 25.4e-3, 1.5)
        segs = segments(only(edge_plots(rendered_plots(tl))))
        @test !isempty(segs)
        @test all(abs(norm(p[[1, 3]]) - 12.7e-3) < 1e-6 for seg in segs for p in seg)

        # edges = false
        for obj in (box, cyl, CubeBeamsplitter(20e-3, λ -> 1.5), System([box, cyl]))
            plots = rendered_plots(obj; edges = false)
            @test isempty(edge_plots(plots))
            @test all(p -> p isa Makie.Mesh, plots)
        end

        # MultiShape: one edge plot per object, after the meshes
        plots = rendered_plots(CubeBeamsplitter(20e-3, λ -> 1.5))
        @test length(edge_plots(plots)) == 1
        @test plots[end] isa Makie.Lines

        # the edges move with the object
        ax = LScene(Figure()[1, 1])
        h = live_render!(ax, box)
        @test h.plots[1] isa Makie.Mesh
        edges = only(edge_plots(h.plots))
        raw = [Vector{Float64}(p) for p in edges[1][] if !any(isnan, p)]
        P0, R0 = h.P0, h.R0
        translate3d!(box, [0.01, -0.02, 0.03])
        yrotate3d!(box, deg2rad(25))
        update_render!(h)
        P, R = BMO.position(box), BMO.orientation(box)
        for i in 1:3
            e = zeros(3)
            e[i] = 1
            @test apply_model(edges, P0 + R0 * e) ≈ P + R * e atol = 1e-6
        end
        @test all(isapprox(apply_model(edges, p), R * R0' * (p - P0) + P; atol = 1e-6) for p in raw)
        remove_render!(h)
    end

    @testset "cemented interfaces" begin
        dl = SphericalDoubletLens(87.9e-3, -105.6e-3, -1000, 6e-3, 3e-3, 25.4e-3, 1.5, 1.6)
        zrotate3d!(dl, deg2rad(35))
        translate3d!(dl, [0.01, -0.02, 0.015])
        plots = rendered_plots(dl)
        @test plots[1] isa Makie.Mesh
        outer, interface = mesh_plots(plots)
        check_material(interface, mats[:interface])
        @test interface.backlight[] == 1
        check_material(outer, mats[:refractive])
        s1, s2 = BMO.shape(dl.front), BMO.shape(dl.back)
        pts = points(interface)
        size = mesh_size(points(outer))
        tol = 1e-9 + 1e-6 * size
        @test !isempty(pts)
        @test maximum(abs(BMO.sdf(s1, p)) for p in pts) < tol
        @test maximum(abs(BMO.sdf(s2, p)) for p in pts) < tol
        # the interface covers the cemented spherical cap once
        area = sum(faces(interface)) do f
            a, b, c = (pts[i] for i in f)
            norm(cross(b - a, c - a)) / 2
        end
        R, r = 105.6e-3, 25.4e-3 / 2
        @test area ≈ 2π * R * (R - sqrt(R^2 - r^2)) rtol = 1e-2
        @test is_closed(points(outer), faces(outer), tol)
        @test length(edge_plots(plots)) == 1
        # an explicit color does not change the interface
        _, interface = mesh_plots(rendered_plots(dl; color = :red))
        check_material(interface, mats[:interface])
        @test isempty(edge_plots(rendered_plots(dl; edges = false)))

        # triplet: both interfaces in one mesh
        tl = SphericalTripletLens(60e-3, 25e-3, -25e-3, -60e-3, 3e-3, 8e-3, 3e-3, 25.4e-3, 1.5, 1.6, 1.5)
        outer, interface = mesh_plots(rendered_plots(tl))
        s1, s2, s3 = BMO.shape(tl.front), BMO.shape(tl.middle), BMO.shape(tl.back)
        pts = points(interface)
        tol = 1e-9 + 1e-6 * mesh_size(points(outer))
        on_12 = [abs(BMO.sdf(s1, p)) < tol && abs(BMO.sdf(s2, p)) < tol for p in pts]
        on_23 = [abs(BMO.sdf(s2, p)) < tol && abs(BMO.sdf(s3, p)) < tol for p in pts]
        @test all(on_12 .| on_23)
        @test any(on_12) && any(on_23)
        @test is_closed(points(outer), faces(outer), tol)

        # other composite objects have no interface mesh
        @test length(mesh_plots(rendered_plots(CubeBeamsplitter(20e-3, λ -> 1.5)))) == 2
        @test length(mesh_plots(rendered_plots(SphericalLens(34.9e-3, -34.9e-3, 6.8e-3, 25.4e-3)))) == 1
    end

    @testset "lighting" begin
        ax = LScene(Figure()[1, 1])
        default_lights = copy(Makie.get_lights(ax.scene))
        # no GLMakie screen: the reduced rig
        @test !Ext._multi_light_backend()
        @test isnothing(studio_lighting!(ax))
        lights = Makie.get_lights(ax.scene)
        @test length(lights) == 1
        @test only(lights) isa Makie.DirectionalLight
        @test only(lights).camera_relative
        @test ax.scene.compute[:ambient_color][] ≈ RGBf(0.45, 0.45, 0.45)
        @test only(lights).color ≈ RGBf(0.75, 0.75, 0.75)
        @test_throws ArgumentError studio_lighting!(ax; preset = :disco)

        # the full rig: ambient + 3 directional lights
        ax = LScene(Figure()[1, 1])
        Ext._apply_lighting!(ax.scene, :studio, true)
        lights = Makie.get_lights(ax.scene)
        @test length(lights) == 3
        @test all(l -> l isa Makie.DirectionalLight && l.camera_relative, lights)
        @test [l.color.r for l in lights] ≈ [0.8, 0.35, 0.3]
        @test ax.scene.compute[:ambient_color][] ≈ RGBf(0.35, 0.35, 0.35)
        @test ax.scene.compute[:lighting_mode][] == Makie.MultiLightShading

        # :none leaves the lights untouched
        ax = LScene(Figure()[1, 1])
        studio_lighting!(ax; preset = :none)
        @test Makie.get_lights(ax.scene) == default_lights

        # live_view
        m = RoundPlanoMirror(25e-3, 5e-3)
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
        gui = live_view(System([m]), beam; detectors = [], lighting = :none)
        @test Makie.get_lights(gui.ax.scene) == default_lights
        close(gui)
        gui = live_view(System([m]), beam; detectors = [])
        @test only(Makie.get_lights(gui.ax.scene)).color ≈ RGBf(0.75, 0.75, 0.75)
        # edges on by default, off via `edges = false`
        @test length(edge_plots(only(gui.system_handles[1].handles).plots)) == 1
        close(gui)
        gui = live_view(System([m]), beam; detectors = [], edges = false)
        @test isempty(edge_plots(only(gui.system_handles[1].handles).plots))
        close(gui)
    end
end

end # module
