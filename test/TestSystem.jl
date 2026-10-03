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

    @testset "Testing push!, pop! and delete!" begin
        m1 = RoundPlanoMirror(25e-3, 5e-3)
        m2 = RoundPlanoMirror(25e-3, 5e-3)
        m3 = RoundPlanoMirror(25e-3, 5e-3)
        m4 = RoundPlanoMirror(25e-3, 5e-3)
        translate3d!(m1, [0, 0.1, 0])
        translate3d!(m2, [0.2, 0, 0])
        translate3d!(m3, [0.3, 0, 0])
        translate3d!(m4, [0.4, 0, 0])
        hit(beam) = BMO.intersection(first(BMO.rays(beam)))
        function solve!(system, beam)
            empty!(beam)
            solve_system!(system, beam)
            return beam
        end

        # the system holds a copy of the vector of objects
        objs = BMO.AbstractObject[m1]
        copied = System(objs)
        push!(copied, m4)
        @test objs == [m1]
        @test System(objs).objects !== objs

        # an empty system exposes no objects and is solved without a hit
        system = System()
        @test isempty(system.objects)
        @test isempty(collect(BMO.objects(system)))
        beam = Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6)
        solve_system!(system, beam)
        @test length(BMO.rays(beam)) == 1
        @test isnothing(hit(beam))

        # push!
        @test push!(system, m1) === system
        @test BMO.object(hit(solve!(system, beam))) === m1
        @test_throws ArgumentError push!(system, m1)
        group = ObjectGroup([m2, m3])
        push!(system, group)
        @test length(system.objects) == 2
        @test system.objects[2] === group
        @test collect(BMO.objects(system)) == [m1, m2, m3]
        # objects of a group that is part of the system
        @test_throws ArgumentError push!(system, m2)
        @test_throws ArgumentError push!(system, ObjectGroup([m3]))
        # nothing is added if one of the objects is invalid
        @test_throws ArgumentError push!(system, m4, m1)
        @test_throws ArgumentError push!(system, m4, m4)
        @test length(system.objects) == 2
        # several objects at once
        other = System()
        push!(other, m1, m4)
        @test collect(BMO.objects(other)) == [m1, m4]

        # delete!
        @test_throws "ObjectGroup" delete!(system, m2)
        @test_throws ArgumentError delete!(system, m2)
        # nothing happens for an object that is not part of the system
        @test delete!(system, m4) === system
        @test length(system.objects) == 2
        @test delete!(system, m1) === system
        @test collect(BMO.objects(system)) == [m2, m3]
        @test isnothing(hit(solve!(system, beam)))

        # a beam group is solved again from its start like a beam
        source = CollimatedSource([0.0, 0, 0], [0.0, 1, 0], 5e-3, 1e-6; num_rings = 2)
        blocked = System([m1])
        solve_system!(blocked, source)
        @test all(b -> length(BMO.rays(b)) == 2, BMO.beams(source))
        delete!(blocked, m1)
        @test empty!(source) === source
        solve_system!(blocked, source)
        @test all(b -> length(BMO.rays(b)) == 1 && isnothing(hit(b)), BMO.beams(source))

        # pop!
        @test pop!(system) === group
        @test isempty(system.objects)
        @test_throws ArgumentError pop!(system)
        @test isnothing(hit(solve!(system, beam)))

        # popat!: the index counts the top-level entries, a group is one of them
        push!(system, m1, group, m4)
        @test_throws BoundsError popat!(system, 0)
        @test_throws BoundsError popat!(system, 4)
        @test length(system.objects) == 3
        @test popat!(system, 2) === group
        @test collect(BMO.objects(system)) == [m1, m4]
        @test popat!(system, 1) === m1
        @test isnothing(hit(solve!(system, beam)))
        @test popat!(system, 1) === m4
        @test isempty(system.objects)
        # a removed object can be added again
        push!(system, m1)
        @test BMO.object(hit(solve!(system, beam))) === m1
    end
end

# A component that absorbs the incoming beam and emits a new one at `shift` from its own
# position, attached with `relaunch!`
struct Reemitter{T, S <: BMO.AbstractShape{T}} <: BMO.AbstractObject{T}
    shape::S
    shift::Vector{T}
end

function Reemitter(shift)
    shape = BMO.QuadraticFlatMesh(0.05)
    zrotate3d!(shape, π)                # normal along −y, like a Detector
    return Reemitter{Float64, typeof(shape)}(shape, shift)
end

const reemit_λ, reemit_w0 = 1e-6, 1e-3

function BMO.interact3d(::BMO.AbstractSystem, r::Reemitter, beam::Beam{T, R},
        ray::R) where {T <: Real, R <: Ray{T}}
    BMO.relaunch!(beam, [Beam(Ray(position(r) + r.shift, direction(ray), BMO.wavelength(ray)))])
    return nothing
end

function BMO.interact3d(::BMO.AbstractSystem, r::Reemitter, g::GaussianBeamlet, ::Int)
    BMO.relaunch!(g,
        [GaussianBeamlet(position(r) + r.shift, direction(g), reemit_λ, reemit_w0)])
    return nothing
end

@testset "relaunch!: re-emitted beams count their path from their own start" begin
    zr, zd, shift = 0.1, 0.3, [2e-3, 0.05, 1e-3]
    function setup()
        r = Reemitter(shift)
        translate3d!(r, [0, zr, 0])
        det = Detector(0.05)
        translate3d!(det, [0, zd, 0] + shift)
        return r, det, System([r, det])
    end
    path = zd - zr                      # from the point of re-emission to the detector

    @testset "Beam" begin
        r, det, system = setup()
        beam = Beam(Ray([0.0, 0, 0], [0.0, 1, 0], reemit_λ))
        solve_system!(system, beam)
        child = only(BMO.children(beam))
        @test isnothing(BMO.AbstractTrees.parent(child))
        @test position(child) ≈ position(r) + shift
        # the path of the parent (zr) is not counted
        @test only(BMO.hits(det)).opl ≈ path
        @test length(child) ≈ path

        # retrace after moving the component: the child is updated in place
        translate3d!(r, [1e-4, 0, 0])
        empty!(det)
        solve_system!(system, beam)
        @test only(BMO.children(beam)) === child
        @test position(child) ≈ position(r) + shift
        @test only(BMO.hits(det)).opl ≈ path

        # another number of beams replaces the children
        new = [Beam(Ray([0.0, 0, 0], [0.0, 1, 0], reemit_λ)) for _ in 1:2]
        BMO.relaunch!(beam, new)
        @test BMO.children(beam) == new
        @test all(isnothing ∘ BMO.AbstractTrees.parent, new)
        BMO.relaunch!(beam, typeof(beam)[])
        @test isempty(BMO.children(beam))
    end

    @testset "GaussianBeamlet" begin
        r, det, system = setup()
        gauss = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], reemit_λ, reemit_w0)
        solve_system!(system, gauss)
        child = only(BMO.children(gauss))
        @test isnothing(BMO.AbstractTrees.parent(child))

        # the field at the detector is that of the same beamlet traced as a source
        reference = GaussianBeamlet(position(r) + shift, [0.0, 1, 0], reemit_λ, reemit_w0)
        ref_det = Detector(0.05)
        translate3d!(ref_det, [0, zd, 0] + shift)
        solve_system!(System([ref_det]), reference)
        hit, ref_hit = only(BMO.hits(det)), only(BMO.hits(ref_det))
        @test hit.l0 == ref_hit.l0 == 0
        for offset in ([0.0, 0, 0], [0.5e-3, 0, 0], [0, 0, -1e-3])
            p = position(det) + offset
            @test BMO.beamlet_hit_field(hit, p) ≈ BMO.beamlet_hit_field(ref_hit, p) rtol = 1e-12
        end

        # retrace after moving the component: the child is updated in place
        translate3d!(r, [1e-4, 0, 0])
        empty!(det)
        solve_system!(system, gauss)
        @test only(BMO.children(gauss)) === child
        @test position(child) ≈ position(r) + shift
    end
end

end # MODULE