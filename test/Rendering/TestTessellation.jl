module TestTessellation

using BeamletOptics
using Makie
using GeometryBasics
using Test
using LinearAlgebra: normalize, norm, dot, cross, I

const BMO = BeamletOptics

const Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

"""Returns the points, normals and faces of the analytic mesh of `s`."""
function mesh_data(s)
    m = Ext._mesh(s)
    return mesh_data(m)
end

function mesh_data(m::GeometryBasics.Mesh)
    pts = [Vector{Float64}(p) for p in GeometryBasics.coordinates(m)]
    nrm = [Vector{Float64}(n) for n in GeometryBasics.normals(m)]
    fcs = [convert.(Int, Tuple(f)) for f in GeometryBasics.faces(m)]
    return pts, nrm, fcs
end

"""Size of the mesh, i.e. the largest width of its bounding box."""
mesh_size(pts) = maximum(maximum(p[i] for p in pts) - minimum(p[i] for p in pts) for i in 1:3)

"""Welds the points by position within `tol`, returns the representative index of each point."""
function weld(pts, tol)
    cell(p) = Tuple(floor.(Int, p ./ tol))
    grid = Dict{NTuple{3, Int}, Vector{Int}}()
    ids = zeros(Int, length(pts))
    for (i, p) in enumerate(pts)
        c = cell(p)
        for d in Iterators.product(-1:1, -1:1, -1:1), j in get(grid, c .+ d, Int[])
            if norm(pts[j] - p) ≤ tol
                ids[i] = ids[j]
                break
            end
        end
        ids[i] == 0 && (ids[i] = i)
        push!(get!(grid, c, Int[]), i)
    end
    return ids
end

"""Returns true if every edge is shared by exactly two faces, after welding the vertices within `tol`."""
function is_closed(pts, fcs, tol)
    ids = weld(pts, tol)
    edges = Dict{Tuple{Int, Int}, Int}()
    for f in fcs
        w = map(i -> ids[i], f)
        allunique(w) || return false
        for k in 1:3
            a, b = w[k], w[mod1(k + 1, 3)]
            key = minmax(a, b)
            edges[key] = get(edges, key, 0) + 1
        end
    end
    return all(==(2), values(edges))
end

centroid(pts, f) = sum(pts[i] for i in f) / 3
face_normal(pts, f) = normalize(cross(pts[f[2]] - pts[f[1]], pts[f[3]] - pts[f[1]]))

"""Checks the primitive `s`: vertices on the surface, orientation, normals and closedness."""
function check_primitive(s)
    pts, nrm, fcs = mesh_data(s)
    size = mesh_size(pts)
    tol = 1e-9 + 1e-6 * size
    @test !isempty(fcs)
    # vertices on the surface
    @test maximum(abs(BMO.sdf(s, p)) for p in pts) < tol
    # outward orientation
    @test all(dot(face_normal(pts, f), sum(nrm[i] for i in f)) > 0 for f in fcs)
    # normals: flat faces at the centroid, curved faces near the vertices (inside the face)
    flat_err, curved_err = 0.0, 0.0
    for f in fcs
        c = centroid(pts, f)
        if all(nrm[i] ≈ nrm[f[1]] for i in f)
            n = BMO.normal3d(s, c)
            flat_err = max(flat_err, norm(n - nrm[f[1]]))
        else
            for i in f
                n = BMO.normal3d(s, pts[i] + 1e-3 * (c - pts[i]))
                curved_err = max(curved_err, acos(clamp(dot(normalize(n), nrm[i]), -1, 1)))
            end
        end
    end
    @test flat_err < 1e-6 # normals are stored as Float32
    @test curved_err < 1e-3
    @test is_closed(pts, fcs, tol)
end

"""
Compares the bounding box of the mesh of `s` with `bounding_box(s)`. The latter probes the SDF
from ±1000 m along the axes, which is only exact if the SDF is exact outside and the extremal
points lie on these axes, hence this check uses the reference pose at the origin. Otherwise the
probe underestimates the extent by `ρ² / 2000`, with `ρ` the distance of the extremal points
from the axis (e.g. the inner radius of a ring).
"""
function check_bounding_box(s)
    pts, _, _ = mesh_data(s)
    tol = 1e-9 + 1e-6 * mesh_size(pts)
    bb = BMO.bounding_box(s)
    for i in 1:3, (sgn, k) in ((-1, 2i - 1), (1, 2i))
        ext = maximum(sgn * p[i] for p in pts)
        ρ = minimum(norm(p[[j for j in 1:3 if j != i]]) for p in pts if sgn * p[i] ≥ ext - tol)
        @test -tol ≤ ext - sgn * bb[k] ≤ tol + ρ^2 / 2000
    end
end

"""
Checks the mesh of the solid with the SDF `solid`: closed and no interior faces, i.e. the outside
of every face is outside of the solid. The offset `ε` exceeds the chord error of the curved faces
(e.g. of the central fan of a lens surface).
"""
function check_solid(pts, nrm, fcs, solid)
    size = mesh_size(pts)
    tol = 1e-9 + 1e-6 * size
    ε = 1e-2 * size
    @test maximum(abs(solid(p)) for p in pts) < tol
    @test all(solid(centroid(pts, f) + ε * face_normal(pts, f)) > 0 for f in fcs)
    @test is_closed(pts, fcs, tol)
end

function posed!(s)
    BMO.zrotate3d!(s, deg2rad(35))
    BMO.xrotate3d!(s, deg2rad(-20))
    BMO.translate3d!(s, [0.01, -0.02, 0.015])
    return s
end

"""Returns the new plots of `render!(ax, x)`."""
function rendered_plots(x)
    ax = LScene(Figure()[1, 1])
    n0 = length(ax.scene.plots)
    render!(ax, x)
    return ax.scene.plots[(n0 + 1):end]
end

# SDF without an analytic mesh, rendered by the marching cubes fallback
mutable struct BlobSDF <: BMO.AbstractSDF{Float64}
    dir::Matrix{Float64}
    transposed_dir::Matrix{Float64}
    pos::Point3{Float64}
end
BlobSDF() = BlobSDF(Matrix(1.0I, 3, 3), Matrix(1.0I, 3, 3), Point3(0.0))
BMO.sdf(s::BlobSDF, p) = norm(BMO._world_to_sdf(s, p)) - 5e-3

@testset "Tessellation" begin
    @test !isnothing(Ext)

    @testset "primitives" begin
        for s in (BMO.BoxSDF(20e-3, 10e-3, 30e-3), BMO.CylinderSDF(12e-3, 4e-3),
                BMO.PlanoSurfaceSDF(5e-3, 25e-3), BMO.RingSDF(10e-3, 4e-3, 6e-3),
                BMO.RightAnglePrismSDF(20e-3, 15e-3), BMO.CutSphereSDF(10e-3, 4e-3),
                BMO.CutSphereSDF(10e-3, -3e-3))
            @testset "$(nameof(typeof(s)))" begin
                check_bounding_box(s)
                check_primitive(posed!(s))
            end
        end
        @testset "SphereSDF" begin
            s = BMO.SphereSDF(8e-3)
            check_bounding_box(s)
            BMO.translate3d!(s, [0.01, -0.02, 0.015])
            check_primitive(s)
        end
        # minimal number of triangles
        @test length(mesh_data(BMO.BoxSDF(1.0, 2.0, 3.0))[3]) == 12
        @test length(mesh_data(BMO.RightAnglePrismSDF(1.0, 2.0))[3]) == 8
    end

    @testset "lenses" begin
        solid(s) = p -> BMO.sdf(s, p)
        lenses = [
            "bi-convex" => SphericalLens(34.9e-3, -34.9e-3, 6.8e-3, 25.4e-3),
            "bi-concave" => SphericalLens(-50e-3, 50e-3, 3e-3, 25.4e-3),
            "plano-convex" => Lens(SphericalSurface(25.8e-3, 25.4e-3), CircularFlatSurface(25.4e-3), 5.3e-3, λ -> 1.5),
            "thin" => ThinLens(34.9e-3, 34.9e-3, 25.4e-3, 1.5),
            "meniscus (left)" => SphericalLens(20e-3, 30e-3, 2e-3, 25.4e-3),
            "meniscus (right)" => SphericalLens(-30e-3, -20e-3, 2e-3, 25.4e-3),
            "bi-concave with ring" => Lens(BMO.BiConcaveLensSDF(50e-3, 50e-3, 3e-3, 20e-3, 25.4e-3), λ -> 1.5),
            "meniscus with ring" => Lens(BMO.MeniscusLensSDF(20e-3, 30e-3, 2e-3, 20e-3, 25.4e-3), λ -> 1.5),
        ]
        for (name, lens) in lenses
            @testset "$name" begin
                posed!(lens)
                plots = rendered_plots(lens)
                @test length(plots) == 1
                @test plots[1] isa Makie.Mesh
                check_solid(mesh_data(plots[1][1][])..., solid(BMO.shape(lens)))
            end
        end
        @test BMO.shape(lenses[5][2]) isa BMO.MeniscusLensSDF
        @test BMO.shape(lenses[6][2]) isa BMO.MeniscusLensSDF

        @testset "doublet" begin
            dl = SphericalDoubletLens(87.9e-3, -105.6e-3, -1000, 6e-3, 3e-3, 25.4e-3, 1.5, 1.6)
            posed!(dl)
            plots = rendered_plots(dl)
            @test length(plots) == 1
            s1, s2 = BMO.shape(dl.front), BMO.shape(dl.back)
            # the cemented interface is interior and removed
            check_solid(mesh_data(plots[1][1][])..., p -> min(BMO.sdf(s1, p), BMO.sdf(s2, p)))
        end
    end

    @testset "other SDF shapes" begin
        # the solids of the cylindric and aspheric lens surfaces
        acyl = Lens(AcylindricalSurface(-15.538e-3, 25e-3, 50e-3, -1.0,
                [0, 1.1926075e-5 * (1e3)^3, -2.9323497e-9 * (1e3)^5]), 7.5e-3, λ -> 1.5)
        for lens in (Lens(CylindricalSurface(5.2e-3, 10e-3, 20e-3), 5.9e-3, λ -> 1.5),
                Lens(CylindricalSurface(-20e-3, 10e-3, 20e-3), CylindricalSurface(20e-3, 10e-3, 20e-3), 3e-3, λ -> 1.5),
                acyl,
                Lens(EvenAsphericalSurface(20e-3, 10e-3, -1.0, [0, 1e3]), 4e-3, λ -> 1.5))
            s = BMO.shape(posed!(lens))
            @test Ext._has_mesh(s)
            pts, _, fcs = mesh_data(s)
            size = mesh_size(pts)
            @test maximum(abs(BMO.sdf(s, p)) for p in pts) < 1e-6 * size
            @test all(BMO.sdf(s, centroid(pts, f) + 1e-3 * size * face_normal(pts, f)) > 0 for f in fcs)
        end
        # conic mirror: front face in the given color, substrate grey
        cm = ConicMirror(200e-3, -0.5, 50e-3; thickness = 6e-3)
        plots = rendered_plots(cm)
        @test length(plots) == 1
        m = plots[1][1][]
        colors = plots[1].color[]
        @test colors isa AbstractVector
        @test length(colors) == length(GeometryBasics.coordinates(m))
        pts, _, fcs = mesh_data(m)
        @test is_closed(pts, fcs, 1e-9 + 1e-6 * mesh_size(pts))
        # mirror with a central hole: DifferenceSDF
        sm = SphericalMirror(200e-3, 6e-3, 25.4e-3; hole_diameter = 5e-3)
        s = BMO.shape(sm)
        @test s isa BMO.DifferenceSDF
        pts, _, fcs = mesh_data(s)
        @test maximum(abs(BMO.sdf(s, p)) for p in pts) < 1e-5
        @test all(norm(centroid(pts, f)[[1, 3]]) > 2.5e-3 * (1 - 1e-3) for f in fcs)
        @test any(norm(centroid(pts, f)[[1, 3]]) < 2.6e-3 for f in fcs) # bore wall
    end

    @testset "objects" begin
        # one mesh plot per color: prisms merged, coating separate
        cbs = CubeBeamsplitter(20e-3, λ -> 1.5)
        plots = rendered_plots(cbs)
        @test length(plots) == 2
        @test all(p -> p isa Makie.Mesh, plots)
        prisms = argmax(p -> length(GeometryBasics.faces(p[1][])), plots)
        pts, _, fcs = mesh_data(prisms[1][])
        @test length(fcs) == 12 # the hypotenuse faces are interior
        @test is_closed(pts, fcs, 1e-9)
        # mesh shapes are flat shaded and lit from both sides
        rr = Retroreflector(25e-3)
        plots = rendered_plots(rr)
        @test length(plots) == 1
        @test plots[1].backlight[] == 1
        pts, nrm, fcs = mesh_data(plots[1][1][])
        @test length(fcs) == 3
        @test all(nrm[i] ≈ face_normal(pts, f) for f in fcs for i in f)
    end

    @testset "marching cubes fallback" begin
        s = BlobSDF()
        @test !Ext._has_mesh(s)
        plots = rendered_plots(s)
        @test length(plots) == 1
        pts, _, _ = mesh_data(plots[1][1][])
        @test !isempty(pts)
        @test maximum(abs(BMO.sdf(s, p)) for p in pts) < 1e-4
        # resolution kwargs are only used by the fallback
        ax = LScene(Figure()[1, 1])
        n0 = length(ax.scene.plots)
        render!(ax, s; x_resolution = 10, y_resolution = 10, z_resolution = 10)
        render!(ax, BMO.BoxSDF(1.0, 1.0, 1.0); x_resolution = 10)
        @test length(ax.scene.plots) == n0 + 2
    end
end

end # module
