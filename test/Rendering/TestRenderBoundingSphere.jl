module TestRenderBoundingSphere

using BeamletOptics
using Makie
using GeometryBasics
using Test
using LinearAlgebra

const BMO = BeamletOptics

const mm = 1e-3

# three closed great circles of 64 segments, separated by two NaN points
const N_SPHERE_POINTS = 3 * 65 + 2

"""Returns the plots that `f()` adds to the `ax`."""
function new_plots(f, ax)
    n0 = length(ax.scene.plots)
    f()
    return ax.scene.plots[(n0 + 1):end]
end

"""Returns the new plots of `render!(ax, x; kwargs...)` in a new `LScene`."""
function rendered_plots(x; kwargs...)
    ax = LScene(Figure()[1, 1])
    return new_plots(() -> render!(ax, x; kwargs...), ax)
end

# The shapes that the solver tests one by one
shapes(obj) = shapes(BMO.shape_trait_of(obj), obj)
shapes(::BMO.SingleShape, obj) = Any[BMO.shape(obj)]
shapes(::BMO.MultiShape, obj) = reduce(vcat, (shapes(part) for part in BMO.shape(obj)); init = Any[])
shapes(group::BMO.AbstractObjectGroup) = reduce(vcat, (shapes(part) for part in BMO.objects(group)); init = Any[])
shapes(sys::System) = reduce(vcat, (shapes(obj) for obj in sys.objects); init = Any[])

# The sphere of a shape. While shape types may still return a tuple (center_local, radius), convert it.
shape_sphere(shape) = shape_sphere(shape, BMO.bounding_sphere_of(shape))
shape_sphere(shape, sphere::BMO.AbstractBoundingSphere) = sphere
shape_sphere(shape, sphere::Tuple) = BMO.SingleBoundingSphere(shape, sphere...)

"""
Returns the bounding spheres of `x` that are drawn, in the order of drawing, derived from
`BMO.bounding_sphere_of`: for a shape its sphere, for an object with a single shape the sphere of the
object, for an object with several parts (also a group) the spheres of the parts and then its own
(main) sphere. Only spheres, no `NoBoundingSphere`.
"""
spheres(x) = filter(s -> !(s isa BMO.NoBoundingSphere), all_spheres(x))
all_spheres(shape::BMO.AbstractShape) = Any[shape_sphere(shape)]
all_spheres(obj::BMO.AbstractObject) = all_spheres(BMO.shape_trait_of(obj), obj)
all_spheres(::BMO.SingleShape, obj) = Any[BMO.bounding_sphere_of(obj)]
function all_spheres(::BMO.MultiShape, obj)
    parts = reduce(vcat, (all_spheres(part) for part in BMO.shape(obj)); init = Any[])
    return push!(parts, BMO.bounding_sphere_of(obj))
end
function all_spheres(group::BMO.AbstractObjectGroup)
    parts = reduce(vcat, (all_spheres(part) for part in BMO.objects(group)); init = Any[])
    return push!(parts, BMO.bounding_sphere_of(group))
end
all_spheres(sys::System) = reduce(vcat, (all_spheres(obj) for obj in sys.objects); init = Any[])

"""Returns the main spheres of `x` that are drawn: those of objects with several parts, also nested ones."""
main_spheres(x) = filter(s -> s isa BMO.MultiBoundingSphere, all_main_spheres(x))
all_main_spheres(::BMO.AbstractShape) = Any[]
all_main_spheres(obj::BMO.AbstractObject) = all_main_spheres(BMO.shape_trait_of(obj), obj)
all_main_spheres(::BMO.SingleShape, obj) = Any[]
function all_main_spheres(::BMO.MultiShape, obj)
    parts = reduce(vcat, (all_main_spheres(part) for part in BMO.shape(obj)); init = Any[])
    return push!(parts, BMO.bounding_sphere_of(obj))
end
function all_main_spheres(group::BMO.AbstractObjectGroup)
    parts = reduce(vcat, (all_main_spheres(part) for part in BMO.objects(group)); init = Any[])
    return push!(parts, BMO.bounding_sphere_of(group))
end
all_main_spheres(sys::System) = reduce(vcat, (all_main_spheres(obj) for obj in sys.objects); init = Any[])

# A shape type without a bounding sphere method
struct NoSphereShape <: BMO.AbstractShape{Float64} end

is_sphere_plot(p) = p isa Makie.Lines && length(p[1][]) == N_SPHERE_POINTS
sphere_plots(plots) = filter(is_sphere_plot, plots)
has_color(p, color) = Makie.to_color(p.color[]) == Makie.to_color(color)
main_plots(plots) = filter(p -> has_color(p, :orange), sphere_plots(plots))

"""Returns the points of the plot `p` in world coordinates, i.e. with its model matrix applied, without the `NaN` separators."""
function world_points(p)
    model = Makie.transformation(p).model[]
    pts = filter(q -> !any(isnan, q), p[1][])
    return map(pts) do q
        w = model * Makie.Point4d(q[1], q[2], q[3], 1.0)
        Makie.Point3d(w[1], w[2], w[3]) ./ w[4]
    end
end

"""Tests that all points of the plot `p` lie on the sphere ``SingleBoundingSphere` in world coordinates."""
function test_on_sphere(p, sphere; rtol = 1e-5)
    center, radius = sphere.pos, sphere.radius
    pts = world_points(p)
    @test length(pts) == 3 * 65
    c = Makie.Point3d(center)
    @test all(q -> isapprox(norm(q - c), radius; rtol), pts)
    return nothing
end

"""
Tests that the plot `p` shows the three great circles of the sphere parallel to
the xy-, yz- and zx-plane, for a plot that has not been moved since it was drawn.
"""
function test_great_circles(p, sphere; rtol = 1e-5)
    test_on_sphere(p, sphere; rtol)
    center, radius = sphere.pos, sphere.radius
    pts = world_points(p)
    c = Makie.Point3d(center)
    for (k, (i, j)) in enumerate(((1, 2), (2, 3), (3, 1)))
        circle = pts[((k - 1) * 65 + 1):(k * 65)]
        n = 6 - i - j
        @test all(q -> isapprox(q[n], c[n]; atol = rtol * radius), circle)
        @test maximum(q -> q[i] - c[i], circle) ≈ radius rtol = rtol
        @test minimum(q -> q[j] - c[j], circle) ≈ -radius rtol = rtol
        @test isapprox(circle[1], circle[end]; atol = rtol * radius)
    end
    return nothing
end

# the plot types, e.g. `mesh` or `lines`
types(plots) = map(Makie.plotfunc, plots)

# A lens whose shape has a bounding sphere on the base branch already, such that the tests do not
# depend on which other shape types have one
meniscus_lens() = Lens(BMO.MeniscusLensSDF(50mm, 80mm, 5mm, 25.4mm), λ -> 1.5)

function test_system()
    lens = SphericalLens(50mm, -50mm, 10mm, 25.4mm, 1.5)
    meniscus = meniscus_lens()
    doublet = SphericalDoubletLens(33.3mm, -22.3mm, -291.1mm, 9mm, 2.5mm, 25.4mm, 1.5, 1.6)
    mirror = RoundPlanoMirror(25.4mm, 5mm)
    cube = CubeBeamsplitter(25mm, λ -> 1.5)
    plate = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, λ -> 1.5)
    filter = RoundPolarizationFilter(10mm)
    polarizer = RoundLinearPolarizer(25.4mm, 1.6mm, 1.6mm, λ -> 1.5)
    prism = RightAnglePrism(20mm, 10mm, 1.5)
    detector = Detector(10mm)
    group = ObjectGroup([RoundPlanoMirror(20mm, 4mm), meniscus_lens()])
    objs = [lens, meniscus, doublet, mirror, cube, plate, filter, polarizer, prism, detector, group]
    for (i, obj) in enumerate(objs)
        translate3d!(obj, [0, i * 50mm, 0])
        zrotate3d!(obj, deg2rad(7 * i))
        xrotate3d!(obj, deg2rad(-3 * i))
    end
    return System(objs)
end

@testset "Bounding sphere rendering" begin
    Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)
    @test !isnothing(Ext)

    @testset "Great circles" begin
        pts = Ext._great_circles(Point3(1.0, 2.0, 3.0), 0.5)
        @test length(pts) == N_SPHERE_POINTS
        @test count(q -> any(isnan, q), pts) == 2
        @test all(isnan, pts[66]) && all(isnan, pts[132])
        @test all(q -> norm(q - Point3(1.0, 2.0, 3.0)) ≈ 0.5, filter(q -> !any(isnan, q), pts))
        @test length(Ext._great_circles(Point3(0.0f0), 1.0f0; segments = 8)) == 3 * 9 + 2
    end

    @testset "Shapes" begin
        ax = LScene(Figure()[1, 1])
        sdf = BMO.MeniscusLensSDF(50mm, 80mm, 5mm, 25.4mm)
        translate3d!(sdf, [10mm, -20mm, 30mm])
        zrotate3d!(sdf, deg2rad(25))
        xrotate3d!(sdf, deg2rad(-40))
        sphere = shape_sphere(sdf)
        @test sphere isa BMO.SingleBoundingSphere

        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, sdf))), ax)
        @test length(plots) == 1
        @test is_sphere_plot(only(plots))
        @test only(plots).linewidth[] == 1
        test_great_circles(only(plots), sphere)

        # keywords: color and linewidth, all others are passed to the plot
        plots = new_plots(ax) do
            BMO.render_bounding_sphere!(ax, sdf; color = :orange, linewidth = 3, visible = false)
        end
        p = only(plots)
        @test p isa Makie.Lines
        @test Makie.to_color(p.color[]) == Makie.to_color(:orange)
        @test p.linewidth[] == 3
        @test p.visible[] == false

        # a shape without a bounding sphere
        none = NoSphereShape()
        @test BMO.bounding_sphere_of(none) isa BMO.NoBoundingSphere
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, none))), ax)
        @test isempty(plots)

        # a mesh: plots follow its bounding sphere
        cube = BMO.CubeMesh(10mm)
        plots = new_plots(() -> BMO.render_bounding_sphere!(ax, cube), ax)
        @test length(plots) == length(spheres(cube))
    end

    @testset "Sphere values" begin
        ax = LScene(Figure()[1, 1])
        sphere = BMO.SingleBoundingSphere(Point3(1.0, 2.0, 3.0), 0.5)
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, sphere))), ax)
        @test length(plots) == 1
        @test is_sphere_plot(only(plots))
        @test has_color(only(plots), :magenta)
        test_great_circles(only(plots), sphere)

        # color and linewidth, `main_color` is not passed to the plot
        plots = new_plots(ax) do
            BMO.render_bounding_sphere!(ax, sphere; color = :red, linewidth = 2, main_color = :green)
        end
        @test only(plots) isa Makie.Lines
        @test has_color(only(plots), :red)
        @test only(plots).linewidth[] == 2

        # no sphere, no plot
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, BMO.NoBoundingSphere()))), ax)
        @test isempty(plots)
    end

    @testset "Objects and systems" begin
        ax = LScene(Figure()[1, 1])

        # SingleShape with a sphere
        lens = meniscus_lens()
        translate3d!(lens, [5mm, 10mm, -15mm])
        yrotate3d!(lens, deg2rad(33))
        @test length(spheres(lens)) == 1
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, lens))), ax)
        test_on_sphere(only(plots), only(spheres(lens)))

        # SingleShape without a sphere
        dummy = NonInteractableObject(joinpath(pkgdir(BeamletOptics), "docs", "src", "assets", "Benchy.stl"))
        @test isempty(spheres(dummy))
        @test isempty(new_plots(() -> BMO.render_bounding_sphere!(ax, dummy), ax))

        # MultiShape, nested: one plot per shape with a sphere, in the order of the shapes. The dummy
        # has no sphere, hence neither the inner nor the outer group has a main sphere.
        group = ObjectGroup([meniscus_lens(), ObjectGroup([meniscus_lens(), dummy]), meniscus_lens()])
        for (i, obj) in enumerate(BMO.Leaves(group))
            translate3d!(obj, [i * 30mm, 0, 0])
        end
        @test length(shapes(group)) == 4
        @test BMO.bounding_sphere_of(group) isa BMO.NoBoundingSphere
        @test length(spheres(group)) == 3
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, group; linewidth = 2))), ax)
        @test length(plots) == length(spheres(group))
        @test isempty(main_plots(plots))
        @test all(p -> p.linewidth[] == 2, plots)
        foreach(test_on_sphere, plots, spheres(group))

        # nested groups with a main sphere each: parts first, the main sphere of a group after its parts
        inner = ObjectGroup([meniscus_lens(), meniscus_lens()])
        translate3d!(inner.objects[2], [0, 30mm, 0])
        outer = ObjectGroup([meniscus_lens(), inner])
        translate3d!(outer.objects[1], [50mm, 0, 0])
        @test BMO.bounding_sphere_of(outer) isa BMO.MultiBoundingSphere
        @test length(main_spheres(outer)) == 2
        plots = new_plots(() -> BMO.render_bounding_sphere!(ax, outer; color = :red, main_color = :green), ax)
        @test length(plots) == length(spheres(outer)) == 5
        @test count(p -> has_color(p, :green), plots) == 2
        @test count(p -> has_color(p, :red), plots) == 3
        foreach(test_on_sphere, plots, spheres(outer))
        @test !has_color(plots[end], :red) # the main sphere of the outer group is drawn last
        # default colors
        plots = new_plots(() -> BMO.render_bounding_sphere!(ax, outer), ax)
        @test length(main_plots(plots)) == 2
        @test count(p -> has_color(p, :magenta), plots) == 3
        # main_color is no plot attribute
        @test !any(p -> haskey(p.attributes, :main_color), plots)

        # system
        sys = test_system()
        expected = spheres(sys)
        @test length(expected) >= 2 # the meniscus lenses
        plots = new_plots(() -> (@test isnothing(BMO.render_bounding_sphere!(ax, sys))), ax)
        @test length(plots) == length(expected)
        @test all(is_sphere_plot, plots)
        foreach(test_great_circles, plots, expected)

        # Axis3
        ax3 = Axis3(Figure()[1, 1])
        plots = new_plots(() -> BMO.render_bounding_sphere!(ax3, sys), ax3)
        @test length(plots) == length(expected)
    end

    @testset "render!(...; show_bounding_sphere)" begin
        sys = test_system()
        expected = spheres(sys)
        @test length(shapes(sys)) > length(sys.objects) # MultiShape objects

        plain = rendered_plots(sys)
        @test isempty(sphere_plots(plain))
        @test types(rendered_plots(sys; show_bounding_sphere = false)) == types(plain)

        plots = rendered_plots(sys; show_bounding_sphere = true)
        found = sphere_plots(plots)
        @test length(found) == length(expected)
        @test types(filter(!is_sphere_plot, plots)) == types(plain)
        # each sphere is drawn once: the plots match the spheres one by one, in the order of the shapes
        foreach(test_on_sphere, found, expected)

        # per object: the same plots as without the keyword, plus one sphere per shape with a sphere
        for obj in sys.objects
            plain = rendered_plots(obj)
            plots = rendered_plots(obj; show_bounding_sphere = true)
            @test types(filter(!is_sphere_plot, plots)) == types(plain)
            @test length(sphere_plots(plots)) == length(spheres(obj))
            @test length(plots) == length(plain) + length(spheres(obj))
            foreach(test_on_sphere, sphere_plots(plots), spheres(obj))
        end

        # other kwargs still reach the plots of the object, not the sphere
        lens = meniscus_lens()
        plots = rendered_plots(lens; show_bounding_sphere = true, color = :red, edges = false)
        @test length(plots) == 2
        @test Makie.to_color(only(filter(p -> p isa Makie.Mesh, plots)).color[]) == Makie.to_color(:red)
        @test length(sphere_plots(plots)) == 1

        # polarizers with their own keywords
        for pol in (RoundPolarizationFilter(10mm), RoundLinearPolarizer(25.4mm, 1.6mm, 1.6mm, λ -> 1.5))
            plain = rendered_plots(pol; show_transmission_axis = false)
            plots = rendered_plots(pol; show_transmission_axis = false, show_bounding_sphere = true)
            @test length(plots) == length(plain) + length(spheres(pol))
            @test length(sphere_plots(plots)) == length(spheres(pol))
        end

        # Axis3
        ax3 = Axis3(Figure()[1, 1])
        plots = new_plots(() -> render!(ax3, sys; show_bounding_sphere = true), ax3)
        @test length(sphere_plots(plots)) == length(expected)
    end

    @testset "live_render!(...; show_bounding_sphere)" begin
        ax = LScene(Figure()[1, 1])

        @testset "SingleShape" begin
            lens = meniscus_lens()
            n_plain = length(rendered_plots(lens))
            h = live_render!(ax, lens; show_bounding_sphere = true)
            @test length(BMO.render_plots(h)) == n_plain + 1
            p = only(sphere_plots(BMO.render_plots(h)))
            test_on_sphere(p, only(spheres(lens)))
            c0 = only(spheres(lens)).pos

            translate3d!(lens, [10mm, -20mm, 30mm])
            update_render!(h)
            @test only(sphere_plots(BMO.render_plots(h))) === p # moved, not drawn again
            @test only(spheres(lens)).pos ≈ c0 + [10mm, -20mm, 30mm]
            test_on_sphere(p, only(spheres(lens)))

            rotate3d!(lens, normalize([1.0, 2.0, -0.5]), deg2rad(50))
            zrotate3d!(lens, deg2rad(-15))
            update_render!(h)
            @test only(sphere_plots(BMO.render_plots(h))) === p
            test_on_sphere(p, only(spheres(lens)))

            n0 = length(ax.scene.plots)
            remove_render!(h)
            @test length(ax.scene.plots) == n0 - n_plain - 1
        end

        @testset "MultiShape" begin
            group = ObjectGroup([meniscus_lens(), meniscus_lens()])
            translate3d!(group.objects[2], [0, 30mm, 0])
            h = live_render!(ax, group; show_bounding_sphere = true)
            found = sphere_plots(BMO.render_plots(h))
            @test length(found) == length(spheres(group)) == 3
            @test length(main_plots(BMO.render_plots(h))) == 1
            foreach(test_on_sphere, found, spheres(group))

            # rigid motion: the plots follow
            translate3d!(group, [5mm, 0, -5mm])
            rotate3d!(group, [0, 0, 1], deg2rad(40))
            update_render!(h)
            @test sphere_plots(BMO.render_plots(h)) == found
            foreach(test_on_sphere, found, spheres(group))

            # a part moved on its own: the object is rendered again, with the spheres
            translate3d!(group.objects[2], [0, 0, 20mm])
            update_render!(h)
            found = sphere_plots(BMO.render_plots(h))
            @test length(found) == length(spheres(group)) == 3
            foreach(test_on_sphere, found, spheres(group))
            remove_render!(h)
        end

        @testset "System" begin
            sys = test_system()
            n0 = length(ax.scene.plots)
            h = live_render!(ax, sys; show_bounding_sphere = true)
            # groups are rendered per member: no main sphere of a group, the system handle has no handle for it
            @test length(sphere_plots(BMO.render_plots(h))) == sum(oh -> length(spheres(BMO.rendered(oh))), BMO.render_children(h))
            for (i, obj) in enumerate(BMO.objects(sys))
                translate3d!(obj, [i * 3mm, 0, -i * 2mm])
                rotate3d!(obj, normalize([1.0, i, 0.3]), deg2rad(11 * i))
            end
            update_render!(h)
            for oh in BMO.render_children(h)
                found = sphere_plots(BMO.render_plots(oh))
                expected = spheres(BMO.rendered(oh))
                @test length(found) == length(expected)
                foreach(test_on_sphere, found, expected)
            end
            remove_render!(h)
            @test length(ax.scene.plots) == n0
        end
    end
end

end
