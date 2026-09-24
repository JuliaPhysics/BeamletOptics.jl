module TestKinematicTrait

using BeamletOptics
using GeometryBasics: Point3
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

"Static test object: opts out of the kinematic API via its trait."
struct Fixed{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end

BMO.kinematic_trait_of(::Fixed) = BMO.Static()

"Static test SDF, used to test the member check of composite SDFs."
mutable struct FixedSDF{T} <: BMO.AbstractSDF{T}
    dir::Matrix{T}
    transposed_dir::Matrix{T}
    pos::Point3{T}
end

FixedSDF() = FixedSDF{Float64}(Matrix(1.0I, 3, 3), Matrix(1.0I, 3, 3), Point3(0.0))

BMO.kinematic_trait_of(::FixedSDF) = BMO.Static()

"Test ray whose beam is static if `fixed` is set, used to test the beam group member check."
mutable struct TaggedRay{T} <: BMO.AbstractRay{T}
    pos::Point3{T}
    dir::Point3{T}
    intersection::BMO.Nullable{BMO.Intersection{T}}
    λ::T
    n::T
    fixed::Bool
end

TaggedRay(fixed::Bool) = TaggedRay{Float64}(Point3(0.0), Point3(0.0, 1.0, 0.0), nothing, 1e-6, 1.0, fixed)

function BMO.kinematic_trait_of(b::BMO.Beam{T, TaggedRay{T}}) where {T}
    return BMO.first_ray(b).fixed ? BMO.Static() : BMO.Movable(BMO.Directed())
end

"Movable test type without kinematic primitives."
struct Bare end

BMO.kinematic_trait_of(::Bare) = BMO.Movable(BMO.Oriented())

"""
Astigmatic beamlets that are declared static. Only `Float32` beamlets are affected, which are
not used anywhere else in the test suite.
"""
const FROZEN = IdSet{Any}()

function BMO.kinematic_trait_of(b::BMO.AstigmaticGaussianBeamlet{Float32})
    return b in FROZEN ? BMO.Static() : BMO.Movable(BMO.Directed())
end

function float32_beamlet()
    c = BMO.Beam(BMO.PolarizedRay(Float32[0, 0, 0], Float32[0, 1, 0], 1.0f-6, ComplexF32[1, 0, 0]))
    r() = BMO.Beam(BMO.Ray(Float32[0, 0, 0], Float32[0, 1, 0], 1.0f-6))
    return BMO.AstigmaticGaussianBeamlet(c, r(), r(), r(), r(), r(), r(), r(), r())
end

"Applies every kinematic verb to `x` and tests that each one throws an `ArgumentError`."
function test_all_verbs_throw(x)
    R = BMO.rotate3d([0.0, 0.0, 1.0], 0.3)
    @test_throws ArgumentError translate3d!(x, [1.0, 0, 0])
    @test_throws ArgumentError translate_to3d!(x, [1.0, 0, 0])
    @test_throws ArgumentError rotate3d!(x, R)
    @test_throws ArgumentError rotate3d!(x, [0.0, 0, 1], 0.3)
    @test_throws ArgumentError rotate3d!(x, R, [1.0, 0, 0])
    @test_throws ArgumentError rotate3d!(x, [0.0, 0, 1], 0.3, [1.0, 0, 0])
    @test_throws ArgumentError xrotate3d!(x, 0.3)
    @test_throws ArgumentError yrotate3d!(x, 0.3)
    @test_throws ArgumentError zrotate3d!(x, 0.3)
    @test_throws ArgumentError align3d!(x, [1.0, 0, 0])
    @test_throws ArgumentError reset_translation3d!(x)
    @test_throws ArgumentError reset_rotation3d!(x)
    return nothing
end

mirror() = RoundPlanoMirror(25mm, 5mm)
lens() = SphericalLens(50mm, -50mm, 5mm, 25mm)
fixed() = Fixed(BMO.CylinderSDF(1mm, 2mm))

"Oriented test subjects, moved away from the identity."
function oriented_subjects()
    subjects = Any[
        BMO.CylinderSDF(1mm, 2mm),
        BMO.CubeMesh(1mm),
        BMO.SphereSDF(1mm) + BMO.CylinderSDF(1mm, 2mm),
        mirror(),
        CubeBeamsplitter(10mm, n -> 1.5),
        ObjectGroup([mirror(), lens()]),
        CollimatedSource([0, -50mm, 0], [0, 1, 0], 5mm, 1e-6; num_rings = 2, num_rays = 40),
    ]
    R = BMO.rotate3d(normalize([1.0, 2.0, 3.0]), 0.7)
    for x in subjects
        rotate3d!(x, R)
        translate3d!(x, [1mm, 2mm, 3mm])
    end
    return subjects
end

"Directed test subjects, moved away from their start state."
function directed_subjects()
    subjects = Any[
        Ray([0, 0, 0], [0, 1, 0]),
        PolarizedRay([0, 0, 0], [0, 1, 0], 1e-6, [1, 0, 0]),
        Beam([0, 0, 0], [0, 1, 0]),
        GaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 1mm),
    ]
    R = BMO.rotate3d(normalize([1.0, 2.0, 3.0]), 0.7)
    for x in subjects
        rotate3d!(x, R)
        translate3d!(x, [1mm, 2mm, 3mm])
    end
    return subjects
end

@testset "Kinematic trait" begin
    @testset "kinematic_trait_of" begin
        oriented = Any[
            BMO.CylinderSDF(1mm, 2mm),
            BMO.CubeMesh(1mm),
            BMO.SphereSDF(1mm) + BMO.CylinderSDF(1mm, 2mm),
            BMO.SphereSDF(1mm) - BMO.CylinderSDF(1mm, 2mm),
            mirror(),
            CubeBeamsplitter(10mm, n -> 1.5),
            ObjectGroup([mirror(), lens()]),
            CollimatedSource([0, -50mm, 0], [0, 1, 0], 5mm, 1e-6; num_rings = 2, num_rays = 40),
        ]
        for x in oriented
            @test BMO.kinematic_trait_of(x) === BMO.Movable(BMO.Oriented())
        end
        directed = Any[
            Ray([0, 0, 0], [0, 1, 0]),
            PolarizedRay([0, 0, 0], [0, 1, 0], 1e-6, [1, 0, 0]),
            Beam([0, 0, 0], [0, 1, 0]),
            GaussianBeamlet([0, 0, 0], [0, 1, 0], 1e-6, 1mm),
        ]
        for x in directed
            @test BMO.kinematic_trait_of(x) === BMO.Movable(BMO.Directed())
        end
        @test BMO.kinematic_trait_of(System([mirror()])) === BMO.Static()
        @test BMO.kinematic_trait_of(1.0) === BMO.Static()
    end

    @testset "Static object" begin
        f = fixed()
        p0 = position(f)
        O0 = copy(BMO.orientation(f))
        test_all_verbs_throw(f)
        @test_throws ArgumentError BMO.direction(f)
        @test position(f) == p0
        @test BMO.orientation(f) == O0
        @test BMO.kinematic_trait_of(f) === BMO.Static()
    end

    @testset "System is static" begin
        system = System([mirror(), lens()])
        test_all_verbs_throw(system)
    end

    @testset "reset_rotation3d! of rays and beams" begin
        @test_throws ArgumentError reset_rotation3d!(Ray([0, 0, 0], [0, 1, 0]))
        @test_throws ArgumentError reset_rotation3d!(Beam([0, 0, 0], [0, 1, 0]))
    end

    @testset "direction and align3d!" begin
        d = [1.0, 2.0, -0.5]
        for x in oriented_subjects()
            @test BMO.direction(x) == BMO.orientation(x)[:, 2]
            align3d!(x, d)
            @test norm(BMO.direction(x) - normalize(d)) < 1e-12
            @test BMO.direction(x) == BMO.orientation(x)[:, 2]
        end
        for x in directed_subjects()
            align3d!(x, d)
            @test norm(BMO.direction(x) - normalize(d)) < 1e-12
        end
    end

    @testset "reset_translation3d! of a traced beam" begin
        m = mirror()
        zrotate3d!(m, deg2rad(45))
        translate3d!(m, [0, 60mm, 0])
        bs = CubeBeamsplitter(10mm, n -> 1.5)
        system = System([bs, m])
        b = Beam([1mm, -50mm, 0.5mm], [0, 1, 0])
        solve_system!(system, b)
        @test length(BMO.rays(b)) > 1
        @test !isempty(BMO.children(b))
        reset_translation3d!(b)
        @test position(b) == zeros(3)
        @test length(BMO.rays(b)) == 1
        @test isnothing(BMO.intersection(first(BMO.rays(b))))
        @test isempty(BMO.children(b))
    end

    @testset "Container member check" begin
        # Object groups, incl. a static group nested in a movable group
        @test_throws ArgumentError ObjectGroup([fixed(), mirror()])
        @test_throws ArgumentError ObjectGroup([mirror(), fixed()])
        @test_throws ArgumentError ObjectGroup([ObjectGroup([fixed(), fixed()]), mirror()])
        @test ObjectGroup([ObjectGroup([fixed(), fixed()]), fixed()]) isa ObjectGroup
        @test_throws ArgumentError ObjectGroup((fixed(), mirror()))
        # Composite SDFs
        s = BMO.SphereSDF(1mm)
        @test_throws ArgumentError FixedSDF() + s
        @test_throws ArgumentError s + FixedSDF()
        @test_throws ArgumentError FixedSDF() - s
        @test_throws ArgumentError s - FixedSDF()
        @test_throws ArgumentError (s + BMO.CylinderSDF(1mm, 2mm)) + FixedSDF()
        @test_throws ArgumentError (FixedSDF() + FixedSDF()) + s
        u = FixedSDF() + FixedSDF()
        @test BMO.kinematic_trait_of(u) === BMO.Static()
        test_all_verbs_throw(u)
        @test BMO.kinematic_trait_of(FixedSDF() - FixedSDF()) === BMO.Static()
        # Beam groups
        mixed = [Beam(TaggedRay(false)), Beam(TaggedRay(true))]
        M = Matrix(1.0I, 3, 3)
        @test_throws ArgumentError PointSource(mixed, 0.1, [0, 0, 0], [0, 1, 0])
        @test_throws ArgumentError PointSource(mixed, 0.1, [0, 0, 0], M)
        @test_throws ArgumentError CollimatedSource(mixed, 1mm, [0, 0, 0], [0, 1, 0])
        @test_throws ArgumentError CollimatedSource(mixed, 1mm, [0, 0, 0], M)
        static_beams = [Beam(TaggedRay(true)), Beam(TaggedRay(true))]
        ps = PointSource(static_beams, 0.1, [0, 0, 0], M)
        @test BMO.kinematic_trait_of(ps) === BMO.Static()
        test_all_verbs_throw(ps)
        movable_beams = [Beam(TaggedRay(false)), Beam(TaggedRay(false))]
        @test BMO.kinematic_trait_of(CollimatedSource(movable_beams, 1mm, [0, 0, 0], M)) ===
              BMO.Movable(BMO.Oriented())
        a_movable = float32_beamlet()
        a_static = float32_beamlet()
        push!(FROZEN, a_static)
        M32 = Matrix(1.0f0I, 3, 3)
        @test_throws ArgumentError AstigmaticBeamGroup([a_movable, a_static], [0, 0, 0], [0, 1, 0])
        @test_throws ArgumentError AstigmaticBeamGroup([a_static, a_movable], [0, 0, 0], M32)
        @test BMO.kinematic_trait_of(AstigmaticBeamGroup([a_movable], [0, 0, 0], M32)) ===
              BMO.Movable(BMO.Oriented())
        # Frames may mix
        @test isnothing(BMO._check_kinematic_members((mirror(), Ray([0, 0, 0], [0, 1, 0]))))
        @test isnothing(BMO._check_kinematic_members(()))
        @test isnothing(BMO._check_kinematic_members([fixed(), System([mirror()])]))
    end

    @testset "Object group trait" begin
        g_static = ObjectGroup([fixed(), fixed()])
        g_movable = ObjectGroup([mirror(), lens()])
        @test BMO.kinematic_trait_of(g_static) === BMO.Static()
        @test BMO.kinematic_trait_of(g_movable) === BMO.Movable(BMO.Oriented())
        @test (@inferred BMO.kinematic_trait_of(g_static)) === BMO.Static()
        @test (@inferred BMO.kinematic_trait_of(g_movable)) === BMO.Movable(BMO.Oriented())
        cs = CollimatedSource([0, -50mm, 0], [0, 1, 0], 5mm, 1e-6; num_rings = 2, num_rays = 40)
        @test (@inferred BMO.kinematic_trait_of(cs)) === BMO.Movable(BMO.Oriented())
        # All verbs throw for a static group, set_pivot3d! still works
        test_all_verbs_throw(g_static)
        set_pivot3d!(g_static, [1.0, 2.0, 3.0])
        @test position(g_static) == [1.0, 2.0, 3.0]
        @test all(position(o) == zeros(3) for o in BMO.objects(g_static))
    end

    @testset "Primitive fallback" begin
        @test_throws ErrorException translate3d!(Bare(), [1.0, 0, 0])
        @test_throws ErrorException rotate3d!(Bare(), Matrix(1.0I, 3, 3))
    end
end

end
