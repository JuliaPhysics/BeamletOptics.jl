module TestStorageShapes

using BeamletOptics
using Test
using Random
using LinearAlgebra: I, norm
const subtypes = BeamletOptics.subtypes

const BMO = BeamletOptics

# Test-only object that holds a shape, so that shapes can go through `save_setup`/`load_setup`
mutable struct ShapeHolder <: BMO.AbstractObject{Float64}
    shape::BMO.AbstractShape{Float64}
end
BMO.to_storage(o::ShapeHolder, ctx) = Dict{String, Any}("shape" => BMO.encode(o.shape, ctx))
BMO.from_storage(::Type{ShapeHolder}, d, ctx) = ShapeHolder(BMO.decode(d["shape"], ctx))
BMO.register_storage_type!(ShapeHolder, "TestShapeHolder")

mutable struct AnyShapeHolder <: BMO.AbstractObject{Float64}
    shape::Any
end
BMO.to_storage(o::AnyShapeHolder, ctx) = Dict{String, Any}("shape" => BMO.encode(o.shape, ctx))
BMO.register_storage_type!(AnyShapeHolder, "TestAnyShapeHolder")

function roundtrip(shapes)
    path = joinpath(mktempdir(), "shapes.bmo")
    BMO.save_setup(path; systems = [System([ShapeHolder(s) for s in shapes])])
    return [o.shape for o in BMO.load_setup(path).systems[1].objects]
end

# Random pose: rotation about a random axis and translation, applied to the whole shape
function place!(rng, s)
    axis = normalize_vec(randn(rng, 3))
    BMO.rotate3d!(s, axis, 2π * rand(rng))
    BMO.translate3d!(s, BMO.Point3((rand(rng, 3) .- 0.5) .* 0.2...))
    return s
end
normalize_vec(v) = v / norm(v)

# Recursively collect all concrete subtypes of `T`
function concrete_subtypes(T)
    out = Any[]
    for S in subtypes(T)
        if isabstracttype(Base.unwrap_unionall(S))
            append!(out, concrete_subtypes(S))
        else
            push!(out, S)
        end
    end
    return out
end


const inch = BMO.inch

function test_shapes()
    [
        ("BoxSDF", BMO.BoxSDF(0.01, 0.02, 0.03)),
        ("CylinderSDF", BMO.CylinderSDF(0.01, 0.02)),
        ("CutSphereSDF", BMO.CutSphereSDF(0.02, 0.005)),
        ("RingSDF", BMO.RingSDF(0.0123, 0.0071, 0.0043)),
        ("RightAnglePrismSDF", BMO.RightAnglePrismSDF(0.02, 0.01)),
        ("PlanoSurfaceSDF", BMO.PlanoSurfaceSDF(0.005, inch)),
        ("SphereSDF", BMO.SphereSDF(0.013)),
        ("ConcaveSphericalSurfaceSDF", BMO.ConcaveSphericalSurfaceSDF(0.05, inch)),
        ("ConvexSphericalSurfaceSDF", BMO.ConvexSphericalSurfaceSDF(0.05, inch)),
        ("MeniscusLensSDF left", BMO.MeniscusLensSDF(0.05, 0.08, 0.005, inch)),
        ("MeniscusLensSDF right", BMO.MeniscusLensSDF(-0.08, -0.05, 0.005, inch)),
        ("MeniscusLensSDF md", BMO.MeniscusLensSDF(0.05, 0.08, 0.005, inch, 1.2inch)),
        ("ConvexAsphericalSurfaceSDF",
            BMO.ConvexAsphericalSurfaceSDF([1.2e-6, -3.4e-9], 0.03, -0.8, 0.02)),
        ("ConcaveAsphericalSurfaceSDF",
            BMO.ConcaveAsphericalSurfaceSDF([1.2e-6, -3.4e-9], -0.03, -0.8, 0.02, 0.025)),
        ("ConvexCylinderSDF", BMO.ConvexCylinderSDF(0.05, 0.02, 0.03)),
        ("ConcaveCylinderSDF", BMO.ConcaveCylinderSDF(0.05, 0.02, 0.03)),
        ("AconvexCylinderSDF", BMO.AconvexCylinderSDF(0.05, 0.02, 0.03, -0.5, [1e-6, 2e-9])),
        ("AconcaveCylinderSDF", BMO.AconcaveCylinderSDF(-0.05, 0.02, 0.03, -0.5, [1e-6, 2e-9])),
        ("OffAxisParaboloidSDF", BMO.OffAxisParaboloidSDF(0.1, 0.03, inch, 0.01)),
        ("UnionSDF ThinLens", BMO.ThinLensSDF(0.05, 0.07)),
        ("UnionSDF BiConvex", BMO.BiConvexLensSDF(0.05, 0.07, 0.008)),
        ("UnionSDF BiConcave md", BMO.BiConcaveLensSDF(0.05, 0.07, 0.004, inch, 1.2inch)),
        ("UnionSDF PlanoConvex", BMO.PlanoConvexLensSDF(0.05, 0.006)),
        ("UnionSDF PlanoConcave md", BMO.PlanoConcaveLensSDF(0.05, 0.004, inch, 1.2inch)),
        ("UnionSDF PlanoConvexAspheric",
            BMO.PlanoConvexAsphericalLensSDF(0.03, 0.008, 0.02, -0.8, [1.2e-6, -3.4e-9])),
        ("UnionSDF PlanoConcaveAspheric",
            BMO.PlanoConcaveAsphericalLensSDF(-0.03, 0.006, 0.02, -0.8, [1.2e-6, -3.4e-9], 0.025)),
        ("UnionSDF nested", BMO.UnionSDF{Float64}(BMO.BiConvexLensSDF(0.05, 0.07, 0.008),
            BMO.BoxSDF(0.01, 0.01, 0.01))),
    ]
end

function cube_mesh()
    V = Float64[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1] .* 0.01
    F = [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; 2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8]
    return BMO.Mesh{Float64}(V, F, BMO.SMatrix{3, 3, Float64, 9}(I), BMO.Point3(0.0, 0, 0), 1.0)
end

function large_mesh(rng)
    n = 600 # 1800 vertex entries, stored as a binary asset
    V = rand(rng, n, 3) .* 0.01
    F = reshape(collect(1:n), :, 3)
    return BMO.Mesh{Float64}(V, F, BMO.SMatrix{3, 3, Float64, 9}(I), BMO.Point3(0.0, 0, 0), 1e-3)
end

@testset "Storage shapes" begin
    @testset "Every SDF type has a tag" begin
        excluded = Any[]  # concrete AbstractSDF subtypes that are intentionally not storable
        types = filter(T -> parentmodule(T) === BeamletOptics,
            vcat(concrete_subtypes(BMO.AbstractSDF), concrete_subtypes(BMO.AbstractMesh)))
        @test !isempty(types)
        for T in types
            T in excluded && continue
            @test (BMO.storage_tag(T); true)
        end
    end

    @testset "SDF round trip" begin
        rng = Xoshiro(20260924)
        cases = test_shapes()
        for (_, s) in cases
            place!(rng, s)
        end
        loaded = roundtrip(last.(cases))
        for ((name, s), l) in zip(cases, loaded)
            @testset "$name" begin
                @test typeof(l) == typeof(s)
                @test maximum(abs, BMO.position(l) - BMO.position(s)) < 1e-14
                @test maximum(abs, BMO.orientation(l) - BMO.orientation(s)) < 1e-14
                c = BMO.position(s)
                Δ = maximum(1:100) do _
                    p = c + BMO.Point3((rand(rng, 3) .- 0.5) .* 0.06...)
                    abs(BMO.sdf(l, p) - BMO.sdf(s, p))
                end
                @test Δ < 1e-12
                if s isa BMO.UnionSDF
                    @test length(l.sdfs) == length(s.sdfs)
                    for (a, b) in zip(l.sdfs, s.sdfs)
                        @test typeof(a) == typeof(b)
                        @test maximum(abs, BMO.position(a) - BMO.position(b)) < 1e-14
                        @test maximum(abs, BMO.orientation(a) - BMO.orientation(b)) < 1e-14
                    end
                end
            end
        end
    end

    @testset "Loaded SDFs keep transposed_dir in sync" begin
        s = place!(Xoshiro(1), BMO.ConvexCylinderSDF(0.05, 0.02, 0.03))
        l = only(roundtrip([s]))
        @test l.transposed_dir == transpose(l.dir)
    end

    @testset "Mesh round trip" begin
        rng = Xoshiro(42)
        small, large = cube_mesh(), large_mesh(rng)
        for m in (small, large)
            BMO.rotate3d!(m, normalize_vec(randn(rng, 3)), 2π * rand(rng))
            BMO.translate3d!(m, BMO.Point3(0.1, -0.2, 0.3))
        end
        path = joinpath(mktempdir(), "mesh.bmo")
        BMO.save_setup(path; systems = [System([ShapeHolder(small), ShapeHolder(large)])])
        ls, ll = [o.shape for o in BMO.load_setup(path).systems[1].objects]
        for (l, m) in ((ls, small), (ll, large))
            @test l isa BMO.Mesh{Float64}
            @test l.vertices == m.vertices
            @test l.faces == m.faces
            @test l.pos == m.pos
            @test l.dir == m.dir
            @test l.scale == m.scale
        end
        # The large mesh is stored as binary assets
        reader = BMO.ZipReader(read(path))
        @test count(n -> startswith(n, "assets/mesh-vertices"), BMO.zip_names(reader)) == 1
        # A loaded mesh can be intersected like the original
        ray = BMO.Ray([0.105, -0.5, 0.305], [0.0, 1.0, 0.0])
        i1, i2 = BMO.intersect3d(small, ray), BMO.intersect3d(ls, ray)
        @test isnothing(i1) == isnothing(i2)
        if !isnothing(i1)
            @test i1.t == i2.t && i1.n == i2.n
        end
    end

    @testset "Only Float64" begin
        path = joinpath(mktempdir(), "f32.bmo")
        @test_throws ArgumentError BMO.save_setup(path;
            systems = [System([AnyShapeHolder(BMO.BoxSDF(1.0f0, 1.0f0, 1.0f0))])])
        @test_throws ArgumentError BMO.save_setup(path;
            systems = [System([AnyShapeHolder(BMO.SphereSDF(1.0f0))])])
    end
end

end # module
