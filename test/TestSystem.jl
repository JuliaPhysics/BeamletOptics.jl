module TestSystem

using BeamletOptics
using Test
# using LinearAlgebra
# using GeometryBasics

const BMO = BeamletOptics

@testset "System" begin
    @testset "Testing implementation" begin
        struct SystemTestBeam{T} <: BMO.AbstractBeam{T, Ray{T}} end
        struct SystemTestObject{T, S} <: BMO.AbstractObject{T} end
        o1 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        o2 = SystemTestObject{Real, BMO.AbstractShape{Real}}()
        system = System(o1)
        beam = SystemTestBeam{Real}()
        # Test missing implementation warnings
        @test_logs (:warn, "Tracing for $(typeof(beam)) not implemented") BMO.trace_system!(
            system,
            beam)
        @test_logs (:warn, "Retracing for $(typeof(beam)) not implemented") BMO.retrace_system!(
            system,
            beam)
    end

    @testset "System metadata and labels" begin
        mutable struct MetadataObject{T} <: BMO.AbstractObject{T}
            id::Int
        end
        object1 = MetadataObject{Float64}(1)
        object2 = MetadataObject{Float64}(2)
        object3 = MetadataObject{Float64}(3)

        system = System(["M1" => object1, "BS" => object2, object3])
        @test_throws ArgumentError System(["M1" => object1, "M1" => object2])
        @test_throws ArgumentError System(["M1" => object1, "M2" => object1])
        BMO.meta(system, object1)[:description] = "front mirror"
        @test BMO.meta(system, object1)[:description] == "front mirror"
        @test BMO.label(system, object1) == "M1"
        @test BMO.label(system, object3) === nothing
        @test system["BS"] === object2

        BMO.label!(system, object1, "M2")
        @test BMO.label(system, object1) == "M2"
        BMO.label!(system, object1, "M2")
        @test BMO.label(system, object1) == "M2"
        @test_throws ArgumentError BMO.label!(system, object2, "M2")
        @test_throws ArgumentError BMO.label(system, MetadataObject{Float64}(4))

        source = "Lens 1"
        BMO.label!(system, object3, SubString(source, 1, 6))
        @test BMO.label(system, object3) == "Lens 1"
        push!(system, "L2" => MetadataObject{Float64}(5))
        @test system["L2"] isa MetadataObject
        @test_throws ArgumentError push!(system, "L2" => MetadataObject{Float64}(6))

        static = StaticSystem(["S1" => object1, object2])
        BMO.meta(static, object1)[:foo] = 42
        @test BMO.meta(static, object1)[:foo] == 42
        @test BMO.label(static, object1) == "S1"
        BMO.label!(static, object2, "S2")
        @test static["S2"] === object2
        copied_static = deepcopy(static)
        @test BMO.meta(copied_static, copied_static["S1"]) !== BMO.meta(static, object1)

        copied = deepcopy(system)
        copied_object = copied["M2"]
        @test BMO.meta(copied, copied_object) !== BMO.meta(system, object1)
        BMO.label!(copied, copied_object, "copy")
        @test BMO.label(system, object1) == "M2"
    end

    # Setup circular multipass cell with flat mirrors
    n_mirrors = 101
    radius = 1
    L = 6 * radius / n_mirrors
    Δθ = 360 / (n_mirrors + 1)
    mirrors = [SquarePlanoMirror2D(L) for _ in 1:n_mirrors]
    θ = 1 * Δθ
    for m in mirrors
        point = radius * [cos(deg2rad(θ)), sin(deg2rad(θ)), 0]
        zrotate3d!(m, deg2rad(θ))
        translate3d!(m, point)
        θ += Δθ
    end
    zrotate3d!.(mirrors, deg2rad(90))

    # Initial ray orientation and position
    dir = [-1, 0, 0]
    Rot = BMO.rotate3d([0, 0, 1], deg2rad(Δθ * 1))
    dir = Vector(Rot * dir)
    origin = [radius, 0, 0] + -1 * dir

    @testset "Testing tracing subroutines" begin
        system = System(mirrors)
        ray = Ray(origin, dir)
        first_obj = mirrors[(n_mirrors + 1) ÷ 2 + 2]
        false_obj = mirrors[(n_mirrors + 1) ÷ 2 + 2 + 1]
        # trace_all
        @test BMO.object(BMO.trace_all(system, ray)) === first_obj
        # trace_one
        @test BMO.object(BMO.trace_one(
            system, ray, BMO.Hint(first_obj))) === first_obj
        @test BMO.object(BMO.trace_one(
            system, ray, BMO.Hint(false_obj))) === first_obj
        # tracing step
        BMO.tracing_step!(system, ray, nothing)
        @test BMO.object(BMO.intersection(ray)) === first_obj
    end

    @testset "Testing system tracing" begin
        system = System(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        # Test trace_system!
        nmax = 10
        BMO.trace_system!(system, beam, r_max = nmax)
        @test length(BMO.rays(beam)) == nmax
        BMO.trace_system!(system, beam, r_max = 1000000)
        @test length(BMO.rays(beam)) == n_mirrors + 1
        first_ray_dir = BMO.direction(first_ray)
        last_ray_dir = BMO.direction(last(BMO.rays(beam)))
        @test 180 - rad2deg(BMO.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test BMO.object(BMO.intersection(first_ray)) ===
              mirrors[(n_mirrors + 1) ÷ 2 + 2]
    end

    @testset "Testing StaticSystem tracing" begin
        # same testset as before
        system = StaticSystem(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        # Test trace_system!
        nmax = 10
        BMO.trace_system!(system, beam, r_max = nmax)
        @test length(BMO.rays(beam)) == nmax
        BMO.trace_system!(system, beam, r_max = 1000000)
        @test length(BMO.rays(beam)) == n_mirrors + 1
        first_ray_dir = BMO.direction(first_ray)
        last_ray_dir = BMO.direction(last(BMO.rays(beam)))
        @test 180 - rad2deg(BMO.angle3d(first_ray_dir, last_ray_dir)) ≈ 2 * Δθ
        @test BMO.object(BMO.intersection(first_ray)) ===
              mirrors[(n_mirrors + 1) ÷ 2 + 2]
    end

    @testset "Testing system retracing" begin
        system = System(mirrors)
        first_ray = Ray(origin, dir)
        beam = Beam(first_ray)
        t1 = @timed BMO.trace_system!(system, beam, r_max = 1000000)
        t2 = @timed BMO.retrace_system!(system, beam) # for precompilation
        t2 = @timed BMO.retrace_system!(system, beam)
        if t1.time < t2.time
            @warn "Retracing took longer than tracing, something might be bugged...\n   Tracing: $(t1.time) s\n   Retracing: $(t2.time) s"
        end
    end
end

end # MODULE