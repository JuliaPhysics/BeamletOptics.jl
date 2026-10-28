module TestSharedParts

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

"MultiShape test object whose parts are given as a tuple, repetitions allowed."
struct Parts{T, P <: Tuple} <: BMO.AbstractObject{T}
    parts::P
end

Parts(parts...) = Parts{Float64, typeof(parts)}(parts)

BMO.shape_trait_of(::Parts) = BMO.MultiShape()
BMO.shape(p::Parts) = p.parts

"SingleShape test object on a given, possibly shared, shape."
struct Single{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end

"Shape of a third party that implements the primitives without type annotations."
struct Foreign{T} <: BMO.AbstractShape{T} end

BMO.translate3d!(::BMO.Movable, ::Foreign, offset) = nothing
BMO.rotate3d!(::BMO.Movable, ::Foreign, R::AbstractMatrix) = nothing

cyl() = BMO.CylinderSDF(1mm, 2mm)

"Returns `x` after moving it to `p`."
placed(x, p) = (translate3d!(x, p); x)

Rz(θ) = BMO.rotate3d([0, 0, 1], θ)

@testset "Shared parts" begin
    @testset "Carried methods and the component interface" begin
        # the carried methods do not collide with a component that implements the primitives
        # without type annotations
        amb = Test.detect_ambiguities(BeamletOptics)
        @test !any(pair -> any(m -> m.module === @__MODULE__, pair), amb)
        @test !any(pair -> any(m -> m.name in (:translate3d!, :rotate3d!), pair), amb)
    end

    @testset "Shape listed twice" begin
        s1, s2 = cyl(), placed(cyl(), [1.0, 0, 0])
        o = Parts(s1, s2, s2)
        translate3d!(o, [0, 0, 1.0])
        @test isapprox(position(s2), [1, 0, 1]; atol = 1e-12)
        # about position(o) = position(s1)
        rotate3d!(o, Rz(π / 2))
        @test isapprox(position(s2), [0, 1, 1]; atol = 1e-12)
        @test isapprox(orientation(s2), Rz(π / 2); atol = 1e-12)
    end

    @testset "Object on a listed shape" begin
        for object_first in (true, false)
            s1, s2 = cyl(), placed(cyl(), [1.0, 0, 0])
            c = Single(s1)
            o = object_first ? Parts(c, s1, s2) : Parts(s1, c, s2)
            translate3d!(o, [0, 0, 1.0])
            @test isapprox(position(s1), [0, 0, 1]; atol = 1e-12)
            @test isapprox(position(s2), [1, 0, 1]; atol = 1e-12)
            rotate3d!(o, Rz(π / 2))
            @test isapprox(position(s1), [0, 0, 1]; atol = 1e-12)
            @test isapprox(orientation(s1), Rz(π / 2); atol = 1e-12)
            @test isapprox(position(s2), [0, 1, 1]; atol = 1e-12)
        end
    end

    @testset "Object in two subgroups" begin
        x = Single(placed(cyl(), [1.0, 0, 0]))
        g = ObjectGroup([ObjectGroup([x]), ObjectGroup([x])])
        translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(x), [1, 0, 1]; atol = 1e-12)
        @test isapprox(position(g), [0, 0, 1]; atol = 1e-12)
        # about the group center, which moved to (0, 0, 1)
        rotate3d!(g, Rz(π / 2))
        @test isapprox(position(x), [0, 1, 1]; atol = 1e-12)
        @test isapprox(orientation(x), Rz(π / 2); atol = 1e-12)
        # derived verbs are built from the primitives and move x once as well
        rotate3d!(g, [0, 0, 1], -π / 2, [0, 0, 1])
        @test isapprox(position(x), [1, 0, 1]; atol = 1e-12)
        @test isapprox(orientation(x), I; atol = 1e-12)
    end

    @testset "Object in a group and a subgroup" begin
        x = Single(placed(cyl(), [1.0, 0, 0]))
        g = ObjectGroup([x, ObjectGroup([x])])
        translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(x), [1, 0, 1]; atol = 1e-12)
    end

    @testset "Nested group with shifted center" begin
        x = Single(placed(cyl(), [1.0, 0, 0]))
        inner = ObjectGroup([x])
        set_pivot3d!(inner, [2.0, 0, 0])
        outer = ObjectGroup([inner])
        rotate3d!(outer, Rz(π / 2))
        @test isapprox(position(inner), [0, 2, 0]; atol = 1e-12)
        @test isapprox(orientation(inner), Rz(π / 2); atol = 1e-12)
        @test isapprox(position(x), [0, 1, 0]; atol = 1e-12)
        @test isapprox(orientation(x), Rz(π / 2); atol = 1e-12)
    end

    @testset "Composite SDF" begin
        # an operand listed twice by the composite itself
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        u = BMO.UnionSDF{Float64}(a, b, b)
        translate3d!(u, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 1]; atol = 1e-12)
        rotate3d!(u, Rz(π / 2))
        @test isapprox(position(b), [0, 1, 1]; atol = 1e-12)
        # a composite that is a part itself shares the set with the other parts
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        o = Parts(a + b, b)
        translate3d!(o, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 1]; atol = 1e-12)
        # the composite is the shape of a SingleShape object, the operand is listed on its own
        for object_first in (true, false)
            # a cylinder, since a sphere ignores rotations
            a, b = BMO.SphereSDF(1.0), placed(cyl(), [1.0, 0, 0])
            c = Single(a + b)
            o = object_first ? Parts(c, b) : Parts(b, c)
            translate3d!(o, [0, 0, 1.0])
            @test isapprox(position(a), [0, 0, 1]; atol = 1e-12)
            @test isapprox(position(b), [1, 0, 1]; atol = 1e-12)
            rotate3d!(o, Rz(π / 2), [0, 0, 1.0])
            @test isapprox(position(b), [0, 1, 1]; atol = 1e-12)
            @test isapprox(orientation(b), Rz(π / 2); atol = 1e-12)
        end
        # the same in a group
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        g = ObjectGroup([Single(a + b), Single(b)])
        translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 1]; atol = 1e-12)
    end

    @testset "Detector in two subgroups" begin
        pd = Detector(50mm)
        translate3d!(pd, [0, 100mm, 0])
        solve_system!(System([pd]), Beam([0.0, 0, 0], [0.0, 1, 0]))
        @test BMO.hit_count(pd) == 1
        g = ObjectGroup([ObjectGroup([pd]), ObjectGroup([pd])])
        translate3d!(g, [0, 0, 1mm])
        @test isnothing(BMO.hits(pd))
        @test isapprox(position(pd), [0, 100mm, 1mm]; atol = 1e-12)
        solve_system!(System([pd]), Beam([0.0, 0, 1mm], [0.0, 1, 0]))
        @test BMO.hit_count(pd) == 1
        D = orientation(pd)
        zrotate3d!(g, 0.1)
        @test isnothing(BMO.hits(pd))
        @test isapprox(orientation(pd), Rz(0.1) * D; atol = 1e-12)
    end
end

end # module
