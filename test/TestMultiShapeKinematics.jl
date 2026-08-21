module TestMultiShapeKinematics

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "MultiShape Kinematics Deduplication" begin
    mutable struct MutableShape{T} <: BMO.AbstractShape{T}
        pos::Point3{T}
        dir::Matrix{T}
    end
    MutableShape(p::Point3{T}) where {T} = MutableShape{T}(p, Matrix{T}(I, 3, 3))
    Base.position(s::MutableShape) = s.pos
    BMO.position!(s::MutableShape, p) = (s.pos = Point3(p))
    BMO.orientation(s::MutableShape) = s.dir
    BMO.orientation!(s::MutableShape, d) = (s.dir = Matrix(d))

    struct MockCoatedOptic{T, S1 <: BMO.AbstractShape{T}, S2 <: BMO.AbstractShape{T}, C <: BMO.AbstractCoating{T}} <: BMO.AbstractObject{T}
        front::S1
        rear::S2
        coating::C
        pos::Base.RefValue{Point3{T}}
        dir::Base.RefValue{Matrix{T}}
    end

    BMO.shape_trait_of(::MockCoatedOptic) = BMO.MultiShape()
    BMO.shape(o::MockCoatedOptic) = (o.front, o.rear, o.coating)
    Base.position(o::MockCoatedOptic) = o.pos[]
    BMO.position!(o::MockCoatedOptic, p) = (o.pos[] = Point3(p))
    BMO.orientation(o::MockCoatedOptic) = o.dir[]
    BMO.orientation!(o::MockCoatedOptic, d) = (o.dir[] = Matrix(d))

    p_front = Point3(0.0, 0.0, 10.0)
    p_rear  = Point3(0.0, 0.0, 15.0)

    s_front = MutableShape(p_front)
    s_rear  = MutableShape(p_rear)

    # front_coating shares the exact same mutable shape instance s_front
    front_coating = Coating(s_front, FresnelInterface())

    optic = MockCoatedOptic(
        s_front,
        s_rear,
        front_coating,
        Ref(Point3(0.0, 0.0, 10.0)),
        Ref(Matrix{Float64}(I, 3, 3))
    )

    @testset "unique_shapes extraction" begin
        uniq = BMO.unique_shapes(optic)
        @test length(uniq) == 2
        @test s_front in uniq
        @test s_rear in uniq
    end

    @testset "Deduplicated translation" begin
        offset = Point3(5.0, -2.0, 3.0)
        translate3d!(optic, offset)

        @test position(s_front) ≈ p_front + offset
        @test position(s_rear)  ≈ p_rear + offset
        @test position(front_coating) ≈ p_front + offset
        @test position(optic) ≈ p_front + offset
    end

    @testset "Deduplicated rotation" begin
        θ = π / 2
        axis = Point3(0.0, 1.0, 0.0)

        initial_front_pos = position(s_front)
        initial_rear_pos  = position(s_rear)

        rotate3d!(optic, axis, θ)

        # Normal vector must rotate by 90 degrees exactly (not 180 degrees from double-rotation)
        R = BMO.rotate3d(axis, θ)
        @test orientation(s_front) ≈ R atol=1e-12
        @test orientation(front_coating) ≈ R atol=1e-12
        @test orientation(optic) ≈ R atol=1e-12

        # Front shape is pivot center, so its position remains unchanged
        @test position(s_front) ≈ initial_front_pos atol=1e-12
        
        # Rear shape is swung around pivot
        expected_rear_offset = R * (initial_rear_pos - initial_front_pos)
        @test position(s_rear) ≈ initial_front_pos + expected_rear_offset atol=1e-12
    end

    @testset "DetectorSurface record! interface" begin
        d_hits = []
        d_surf = DetectorSurface(d_hits)
        ray = Ray(Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 1.0))
        hit = Intersection(10.0, Point3(0.0, 0.0, 10.0), Point3(0.0, 0.0, -1.0))
        trans = Transition(Ambient(), Ambient(), true)

        res = BMO.interact3d(d_surf, trans, hit, ray)
        @test res === nothing
        @test length(d_hits) == 1
        @test d_hits[1][1] ≈ Point3(0.0, 0.0, 10.0)
    end
end

end # module
