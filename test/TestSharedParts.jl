module TestSharedParts

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

"MultiShape test object whose shapes are given as a tuple, repetitions allowed."
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

"SingleShape test object that opts out of the kinematic API."
struct Fixed{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
end

BMO.kinematic_trait_of(::Fixed) = BMO.Static()

"SingleShape test object that counts the calls of its primitives, like a detector that empties itself."
mutable struct Counting{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    const shape::S
    calls::Int
end

Counting(s) = Counting(s, 0)

function BMO.translate3d!(::BMO.Movable, c::Counting, offset)
    c.calls += 1
    BMO.translate3d!(BMO.shape_trait_of(c), c, offset)
    return nothing
end

function BMO.rotate3d!(::BMO.Movable, c::Counting, R::AbstractMatrix)
    c.calls += 1
    BMO.rotate3d!(BMO.shape_trait_of(c), c, R)
    return nothing
end

cyl() = BMO.CylinderSDF(1mm, 2mm)

"Returns `x` after moving it to `p`."
placed(x, p) = (translate3d!(x, p); x)

Rz(θ) = BMO.rotate3d([0, 0, 1], θ)

@testset "Shared parts" begin
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

    @testset "Subgroup in two groups" begin
        x = Single(placed(cyl(), [1.0, 0, 0]))
        sub = ObjectGroup([x])
        g = ObjectGroup([ObjectGroup([sub]), sub])
        translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(x), [1, 0, 1]; atol = 1e-12)
        @test isapprox(position(sub), [0, 0, 1]; atol = 1e-12)
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

    @testset "Empty subgroup" begin
        empty = ObjectGroup([])
        set_pivot3d!(empty, [1.0, 0, 0])
        g = ObjectGroup([empty, Single(cyl())])
        rotate3d!(g, Rz(π / 2))
        @test isapprox(position(empty), [0, 1, 0]; atol = 1e-12)
        @test isapprox(orientation(empty), Rz(π / 2); atol = 1e-12)
    end

    @testset "Composite SDF" begin
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        u = BMO.UnionSDF{Float64}(a, b, b)
        translate3d!(u, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 1]; atol = 1e-12)
        rotate3d!(u, Rz(π / 2))
        @test isapprox(position(b), [0, 1, 1]; atol = 1e-12)
        # an operand that is also the shape of another object would move twice
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        g = ObjectGroup([Single(a + b), Single(b)])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 0]; atol = 1e-12)
        # the same within one composite
        a, b = BMO.SphereSDF(1.0), placed(BMO.SphereSDF(0.5), [1.0, 0, 0])
        u = BMO.UnionSDF{Float64}(a + b, b)
        @test_throws ArgumentError translate3d!(u, [0, 0, 1.0])
        @test isapprox(position(b), [1, 0, 0]; atol = 1e-12)
    end

    @testset "Primitives of an object type apply in nested groups" begin
        c = Counting(cyl())
        g = ObjectGroup([ObjectGroup([c]), ObjectGroup([c])])
        translate3d!(g, [0, 0, 1.0])
        @test c.calls == 1
        @test isapprox(position(c), [0, 0, 1]; atol = 1e-12)
        zrotate3d!(g, π / 2)
        # a rotation about the group pivot is derived from rotate3d! and translate3d!
        @test c.calls == 3
        o = Parts(cyl())
        g = ObjectGroup([o, ObjectGroup([c])])
        translate3d!(g, [0, 0, 1.0])
        @test c.calls == 4
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

    @testset "Errors before anything moves" begin
        # a static object as a part of a MultiShape object
        s = cyl()
        o = Parts(s, Fixed(cyl()))
        @test_throws ArgumentError translate3d!(o, [0, 0, 1.0])
        g = ObjectGroup([o])
        @test_throws ArgumentError zrotate3d!(g, π / 2)
        @test isapprox(position(s), [0, 0, 0]; atol = 1e-12)
        @test isapprox(orientation(s), I; atol = 1e-12)
        @test isapprox(orientation(g), I; atol = 1e-12)
        # two objects on the same shape; two `Single(shared)` would be one object, since an
        # immutable struct with the same fields is `===`
        x = Single(placed(cyl(), [1.0, 0, 0]))
        shared = cyl()
        g = ObjectGroup([x, ObjectGroup([Single(shared), Counting(shared)])])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test_throws ArgumentError zrotate3d!(g, π / 2)
        @test isapprox(position(x), [1, 0, 0]; atol = 1e-12)
        @test isapprox(position(shared), [0, 0, 0]; atol = 1e-12)
        @test isapprox(position(g), [0, 0, 0]; atol = 1e-12)
        # a MultiShape object and a SingleShape object on one of its shapes
        s = cyl()
        g = ObjectGroup([Parts(s, cyl()), Single(s)])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(s), [0, 0, 0]; atol = 1e-12)
        # the shape is one level down, below an object part
        s = cyl()
        g = ObjectGroup([Parts(Single(s), cyl()), Single(s)])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(s), [0, 0, 0]; atol = 1e-12)
        # within one MultiShape object
        s = cyl()
        @test_throws ArgumentError translate3d!(Parts(Single(s), s), [0, 0, 1.0])
        @test isapprox(position(s), [0, 0, 0]; atol = 1e-12)
        # a cemented doublet and its own front lens
        dl = SphericalDoubletLens(40mm, -30mm, -120mm, 6mm, 2.5mm, 25mm, _ -> 1.62, _ -> 1.72)
        p_front, p_back = position(dl.front), position(dl.back)
        g = ObjectGroup([dl, dl.front])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(dl.front), p_front; atol = 1e-12)
        @test isapprox(position(dl.back), p_back; atol = 1e-12)
        # a static object two levels down
        a = Single(cyl())
        g = ObjectGroup([a, Parts(cyl(), Parts(cyl(), Fixed(cyl())))])
        @test_throws ArgumentError translate3d!(g, [0, 0, 1.0])
        @test isapprox(position(a), [0, 0, 0]; atol = 1e-12)
        @test isapprox(position(g), [0, 0, 0]; atol = 1e-12)
    end

    @testset "Built-in MultiShape objects in a group" begin
        dl = SphericalDoubletLens(40mm, -30mm, -120mm, 6mm, 2.5mm, 25mm, _ -> 1.62, _ -> 1.72)
        cbs = CubeBeamsplitter(10mm, n -> 1.5)
        pbs = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, n -> 1.5)
        p = [position(x) for x in (dl.front, dl.back, cbs, pbs)]
        g = ObjectGroup([ObjectGroup([dl, cbs]), ObjectGroup([pbs, dl])])
        translate3d!(g, [0, 0, 1.0])
        q = [position(x) for x in (dl.front, dl.back, cbs, pbs)]
        @test all(isapprox.(q .- p, Ref([0, 0, 1.0]); atol = 1e-12))
    end
end

end # module
