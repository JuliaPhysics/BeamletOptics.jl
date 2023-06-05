using SCDI
using Test
using LinearAlgebra
using UUIDs

@testset "Utilities" begin
    @testset "Testing euclidic norm utilites" begin
        v = [1.0, 0.0, 1.0]
        @test isapprox(SCDI.norm3d(v), sqrt(2))
        @test isapprox(SCDI.normalize3d(v), [sqrt(2) / 2, 0, sqrt(2) / 2])
        SCDI.normalize3d!(v)
        @test isapprox(v, [sqrt(2) / 2, 0, sqrt(2) / 2])
    end

    @debug "Testing orthogonal3d for right hand rule and unit length"
    orth = SCDI.orthogonal3d([2, 0, 0], [0, 0, 1])
    @test isapprox(orth, [0, -1, 0])

    @debug "Testing rotate3d for clockwise dir. and conservation of length"
    Rot = SCDI.rotate3d([0, 0, 1], π / 2)
    @test isapprox(Rot * [1, 0, 0], [0, 1, 0])

    @debug "Testing align3d for rotation and conservation of length"
    target = [0, 0, 1]
    Rot = SCDI.align3d([1, 0, 0], target)
    @test isapprox(Rot * [1, 0, 0], target)

    @debug "Testing angle3d for resulting angle"
    a = SCDI.angle3d([1, 0, 0], [0, 0, 1])
    @test isapprox(a, π / 2)

    @debug "Testing fast_dot3d"
    c = SCDI.fast_dot3d([1, 2, 3], [4, 5, 6])
    @test c == 32

    @debug "Testing fast_cross3d"
    c = SCDI.fast_cross3d([1, 0, 0], [0, 1, 0])
    @test isapprox(c, [0, 0, 1])

    @debug "Testing fast_cross3d!"
    c = zeros(Int, 3)
    SCDI.fast_cross3d!(c, [1, 0, 0], [0, 1, 0])
    @test isapprox(c, [0, 0, 1])

    @debug "Testing fast_sub3d!"
    c = zeros(Int, 3)
    SCDI.fast_sub3d!(c, [4, 5, 6], [1, 2, 3])
    @test isapprox(c, [3, 3, 3])

    @debug "Testing line_point_distance3d"
    pos = [0, 0, 0]
    dir = [1, 0, 0]
    point = [5, 1, 1]
    d = SCDI.line_point_distance3d(pos, dir, point)
    @test isapprox(d, √2)

    @testset "Testing reflection3d" begin
        for dx in -1:1, dy in -1:1
            @test isapprox(SCDI.reflection3d([dx, dy, 1], [0, 0, -1]), [dx, dy, -1])
        end
    end

    @testset "Testing refraction3d" begin
        normal = [0, 0, 1]
        @testset "Test from vacuum into medium" begin
            n1 = 1.0
            n2 = 1.5
            for θ1 in 0:π/8:π/2
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out = SCDI.refraction3d(dir_in, normal, n1, n2)
                θ2 = SCDI.angle3d(-normal, dir_out)
                # 2D-equation for refraction validation
                θ3 = asin(n1 / n2 * sin(θ1))
                @test isapprox(θ2, θ3)
            end
        end
        @testset "Test from medium into vacuum" begin
            n1 = 1.5
            n2 = 1.0
            for θ1 in 0:π/8:π/2
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out = SCDI.refraction3d(dir_in, normal, n1, n2)
                if θ1 > asin(n2 / n1)
                    # Test for total reflection
                    θ2 = SCDI.angle3d(dir_out, normal)
                    @test isapprox(θ1, θ2)
                else
                    # Test for refraction
                    θ2 = SCDI.angle3d(-normal, dir_out)
                    θ3 = asin(n1 / n2 * sin(θ1))
                    @test isapprox(θ2, θ3)
                end
            end
        end
    end
end

@testset "Types" begin
    @debug "Testing abstract type definitions"
    @test isdefined(SCDI, :AbstractEntity)
        @test isdefined(SCDI, :AbstractObject)
            @test isdefined(SCDI, :AbstractMirror)
            @test isdefined(SCDI, :AbstractPrism)
        @test isdefined(SCDI, :AbstractShape)
            @test isdefined(SCDI, :AbstractMesh)
        @test isdefined(SCDI, :AbstractRay)

    # Generate test structs
    struct Shapeless{T} <: SCDI.AbstractShape{T}
        pos::Vector{T}
        dir::Matrix{T}
    end

    struct Ray{T} <: SCDI.AbstractRay{T}
        pos::Vector{T}
        dir::Vector{T}
    end

    @testset "Testing abstract type interfaces" begin
        object = Shapeless([0,0,0], Matrix{Int}(I, 3, 3))
        SCDI.translate3d!(object, [1,1,1])
        SCDI.position!(object, [2,2,2])
        SCDI.reset_translation3d!(object)
        @test SCDI.position(object) == object.pos
        @test SCDI.position(object) == zeros(3)
        # The following test are expected to do nothing but not throw exceptions
        SCDI.rotate3d!(object, [0,0,1], π)
        @test all(SCDI.orientation(object)[[1, 5]] .== -1) && SCDI.orientation(object)[9] == 1
        SCDI.xrotate3d!(object, π)
        SCDI.yrotate3d!(object, π)
        SCDI.zrotate3d!(object, π)
        SCDI.reset_rotation3d!(object)
        @test SCDI.orientation(object) == Matrix{Int}(I, 3, 3)
        @test_logs (:warn,) SCDI.intersect3d(object, SCDI.Ray([0,0,1],[0,0,1]))
    end

    @testset "Testing abstract ray" begin
        object = Shapeless([1,0,0], Matrix{Int}(I, 3, 3))
        ray = Ray([0,0,0], [1,0,0])
        @test SCDI.isinfrontof(object, ray) == true
        SCDI.direction!(ray, -[1,0,0])
        @test SCDI.isinfrontof(object, ray) == false
        SCDI.direction!(ray, [0,1,0])
        @test SCDI.isinfrontof(object, ray) == false
        SCDI.direction!(ray, [1,1,0])
        @test SCDI.isinfrontof(object, ray) == true
    end
end

@testset "Rays" begin
    @debug "Testing Ray struct definition"
    @test isdefined(SCDI, :Ray)
    pos = [0, 0, 0]
    dir = [1, 1, 1]
    ray = SCDI.Ray(pos, dir)
    @test ismutable(ray)
    @test eltype(ray.pos) == Float64

    @debug "Testing Ray constructor (dir normalization and init. Inf length)"
    @test isapprox(SCDI.norm3d(ray.dir), 1)
    @test isnothing(ray.intersection)

    @debug "Testing Beam struct definition"
    @test isdefined(SCDI, :Beam)
    beam = SCDI.Beam([ray])
    @test ismutable(beam)
end

@testset "Mesh" begin
    # NOTE: the "Mesh" testset is mutating. Errors/fails might lead to subsequent tests failing too!
    @testset "Testing type definitions" begin
        @debug "Testing Mesh struct definition"
        @test isdefined(SCDI, :Mesh)
    end

    # Generate cube since types are defined
    include("Cube.jl")

    @debug "Generating moving test cube with scale 1"
    foo = Cube(1)

    @testset "Testing AbstractMesh getters" begin
        @test typeof(foo) == SCDI.Mesh{Float64}
        @test SCDI.vertices(foo) == foo.vertices
        @test SCDI.faces(foo) == foo.faces
        @test SCDI.orientation(foo) == foo.dir
        @test SCDI.position(foo) == foo.pos
        @test SCDI.scale(foo) == foo.scale
    end

    @testset "Testing translate3d!" begin
        SCDI.translate3d!(foo, -0.5 * [1, 1, 1]) # move COG to origin
        @test minimum(SCDI.vertices(foo)[:, 1]) == -0.5
        @test minimum(SCDI.vertices(foo)[:, 2]) == -0.5
        @test minimum(SCDI.vertices(foo)[:, 3]) == -0.5
        @test maximum(SCDI.vertices(foo)[:, 1]) == 0.5
        @test maximum(SCDI.vertices(foo)[:, 2]) == 0.5
        @test maximum(SCDI.vertices(foo)[:, 3]) == 0.5
        @test all(SCDI.position(foo) .== -0.5)
    end

    @testset "Testing set_new_origin3d!" begin
        SCDI.set_new_origin3d!(foo)
        @test SCDI.position(foo) == zeros(3)
    end

    @testset "Testing x/y/zrotate3d!" begin
        @testset "Testing xrotate3d!" begin
            SCDI.xrotate3d!(foo, π / 4)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), √2 / 2)
        end

        @testset "Testing yrotate3d!" begin
            SCDI.yrotate3d!(foo, π / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 0.5)
        end

        @testset "Testing zrotate3d!" begin
            SCDI.zrotate3d!(foo, π / 4)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -0.5)
            @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), 0.5)
            @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 0.5)
        end

        @debug "Testing orientation of dir matrix"
        @test all(SCDI.orientation(foo)[[2, 6, 7]] .== 1)
    end

    @testset "Testing reset_rotation3d!" begin
        SCDI.reset_rotation3d!(foo)
        SCDI.reset_translation3d!(foo) # sneaky, useless command for that sweet code coverage
        @test all(SCDI.orientation(foo) .== Matrix{Float64}(I, 3, 3))
    end

    @testset "Testing orthogonal3d" begin
        normal = SCDI.orthogonal3d(foo, 1)
        @test isapprox(normal, [0, 0, -1])
    end

    @testset "Testing scale3d!" begin
        SCDI.scale3d!(foo, 2)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 1]), -1)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 2]), -1)
        @test isapprox(minimum(SCDI.vertices(foo)[:, 3]), -1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 1]), 1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 2]), 1)
        @test isapprox(maximum(SCDI.vertices(foo)[:, 3]), 1)
        @test SCDI.scale(foo) == 2
    end
end

@testset "Intersections" begin
    @testset "Testing Moeller-Trumbore algorithm" begin
        @debug "Testing functor definition"
        @test isdefined(SCDI, :MoellerTrumboreAlgorithm)
        @test hasfield(SCDI.MoellerTrumboreAlgorithm, :kϵ)
        @test hasfield(SCDI.MoellerTrumboreAlgorithm, :lϵ)
        @test isdefined(SCDI, :ray_triangle_intersection)
        @test isconst(SCDI, :ray_triangle_intersection)
        @debug "Defining face in the x-y-plane with z-offset t"
        t = 5
        face = [
            1 1 t
            -1 1 t
            0 -1 t
        ]
        @debug "Defining ray at origin pointing along z-axis"
        pos = [0, 0, 0]
        dir = [0, 0, 1]
        ray = SCDI.Ray(pos, dir)
        # Preallocate memory
        E1 = similar(ray.dir)
        E2 = similar(ray.dir)
        Pv = similar(ray.dir)
        Tv = similar(ray.dir)
        Qv = similar(ray.dir)
        @debug "Testing intersection"
        @test isapprox(SCDI.ray_triangle_intersection(face, ray, E1, E2, Pv, Tv, Qv), t)
        # Check allocations (WARNING: function must have been compiled once for before this test!)
        alloc = @allocated SCDI.ray_triangle_intersection(face, ray, E1, E2, Pv, Tv, Qv)
        if alloc > 16
            @warn "Allocated number of bytes for MTA larger than expected!" alloc
        end
    end

    @testset "Testing intersect3d for meshes" begin
        @debug "Generating test cube"
        t = 5
        s = 1 # scale/2
        foo = Cube(2 * s)
        # Move cube COG to origin
        SCDI.translate3d!(foo, -[s, s, s])
        SCDI.set_new_origin3d!(foo)
        # Align cube edge at t units from origin
        SCDI.translate3d!(foo, [t + s, 0, 0])
        @debug "Defining ray at origin with variable dir"
        pos = [0, 0, 0]
        steps = 10
        for z in -s:(s/steps):s
            # Ray constructed each time for unit-length dir
            dir = [t, 0, z]
            ray = SCDI.Ray(pos, dir)
            @test isapprox(SCDI.intersect3d(foo, ray).t, sqrt(t^2 + z^2))
        end
    end
end

@testset "Interactions" begin
    @test isdefined(SCDI, :Interaction)
    @test isdefined(SCDI, :Mirror)
    @test isdefined(SCDI, :Prism)
end

@testset "System" begin
    # Set up system of four mirrors boxing in the origin in the x-y-plane
    c1 = Cube(2)
    c2 = Cube(2)
    c3 = Cube(2)
    c4 = Cube(2)
    SCDI.translate3d!.([c1, c2, c3, c4], [-[1,1,1], -[1,1,1], -[1,1,1], -[1,1,1]])
    SCDI.translate3d!.([c1, c2, c3, c4], [[2,0,0], -[2,0,0], [0,2,0], -[0,2,0]])
    SCDI.set_new_origin3d!.([c1, c2, c3, c3])
    m1 = SCDI.Mirror(uuid4(), c1)
    m2 = SCDI.Mirror(uuid4(), c2)
    m3 = SCDI.Mirror(uuid4(), c3)
    m4 = SCDI.Mirror(uuid4(), c4)
    system = SCDI.System([m1, m2, m3, m4])
    # Set up ray that will be reflected at 45° angle in the middle of each mirror edge
    ray = SCDI.Ray([0,-1,0], [1,1,0])
    beam = SCDI.Beam([ray])
    # Trace system
    numEl = 8
    alloc1 = @allocated SCDI.solve_system!(system, beam, r_max=8)
    temp = SCDI.Beam(copy(beam.rays))
    # Retrace system
    alloc2 =  @allocated SCDI.solve_system!(system, beam, r_max=8)
    @test SCDI.position.(beam.rays) == SCDI.position.(temp.rays)
    @test length(beam.rays) == 8
    for pos in SCDI.position.(beam.rays)
        @test pos[3] == 0
    end
    if alloc2 >= alloc1
        @warn "No retracing utilized"
    end
    # Remove object from beam path and retrace
    SCDI.translate3d!(c2, [0,0,2])
    SCDI.solve_system!(system, beam, r_max=8)
    @test length(beam.rays) == 3

    #=
    WIP:
    - test for object ID correctness
    - test for retracing correctness
    - remove alloc warning?
    =#
end

@testset "SDFs" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :AbstractSDF)
        @test isdefined(SCDI, :SphereSDF)
        @test isdefined(SCDI, :CylinderSDF)
        @test isdefined(SCDI, :CutSphereSDF)
        @test isdefined(SCDI, :ThinLensSDF)
    end

    @testset "Thin lens focal length" begin
        # define thin lens
        R1 = 1
        R2 = 1
        nl = 1.5
        tl = SCDI.ThinLensSDF(R1, R2, 0.1)
        p = SCDI.Prism(uuid4(), tl, x -> 1.5)
        system = SCDI.System([p])
        
        # compare numerical and analytical focal length
        f_analytical = SCDI.lensmakers_eq(R1, -R2, nl)
        zs = -0.04:0.01:0.04
        for (i, z) in enumerate(zs)
            # skip optical axis ray
            if z ≈ 0
                continue
            end
            xs = 0.1:0.1:1.5
            df = zeros(Float64, length(xs))
            ray = SCDI.Ray([0, -0.5, z], [0, 1, 0], 1e3)
            beam = SCDI.Beam(ray)
            SCDI.solve_system!(system, beam)
            # test if numerical and analytical focal length agree
            for (i, x) in enumerate(xs)
                df[i] = SCDI.line_point_distance3d(beam.rays[end], [0, x, 0])
            end
            @test xs[findmin(df)[2]] ≈ f_analytical
        end
    end
end

@testset "Gaussian beamlet" begin
    @testset "Testing type definitions" begin
        @test isdefined(SCDI, :GaussianBeamlet)
    end

    @testset "Testing analytical equations" begin
        λ = 500e-9
        w0 = 1e-3
        M2 = 1
        zR = SCDI.rayleigh_range(λ, w0, M2)
        # Test Rayleigh range and div. angle against Paschotta
        @test isapprox(zR, 6.28, atol=1e-2)
        @test isapprox(SCDI.beam_waist(zR, w0, zR), sqrt(2)*w0)
        @test isapprox(SCDI.gouy_phase(zR, zR), -π/4)
        @test isapprox(SCDI.wavefront_curvature(zR, zR), 2*zR)
        @test isapprox(SCDI.divergence_angle(λ, w0, M2), 159e-6, atol=1e-6)
    end

    @testset "Testing parameter correctness" begin
        # Gauss beam parameters
        y = -5:0.01:5
        λ_1 = 500
        λ_2 = 1000
        w0_1 = 1
        w0_2 = 2
        M2_1 = 1
        M2_2 = 2
        gauss_1 = SCDI.GaussianBeamlet(SCDI.Ray([0, 0, 0], [0, 1, 0]), λ_1, w0_1, M2 = M2_1)
        gauss_2 = SCDI.GaussianBeamlet(SCDI.Ray([0, 0, 0], [0, 1, 0]), λ_2, w0_2, M2 = M2_2)
        # Calculate analytical values
        zr_1 = SCDI.rayleigh_range(λ_1*1e-9, w0_1*1e-3, M2_1)
        zr_2 = SCDI.rayleigh_range(λ_2*1e-9, w0_2*1e-3, M2_2)
        wa_1 = SCDI.beam_waist.(y, w0_1*1e-3, zr_1)
        wa_2 = SCDI.beam_waist.(y, w0_2*1e-3, zr_2)
        ψa_1 = SCDI.gouy_phase.(y, zr_1)
        ψa_2 = SCDI.gouy_phase.(y, zr_2)
        # Calculate numerical values
        wn_1, Rn_1, ψn_1 = SCDI.gauss_parameters(gauss_1, y)
        wn_2, Rn_1, ψn_2 = SCDI.gauss_parameters(gauss_2, y)
        # Compare beam diameter within 0.1 nm
        @test all(isapprox.(wa_1, wn_1, atol=1e-10))
        @test all(isapprox.(wa_2, wn_2, atol=1e-10))
        # Compare Gouy phase within 0.1 μrad
        @test all(isapprox.(ψa_1, ψn_1, atol=1e-7))
        @test all(isapprox.(ψa_2, ψn_2, atol=1e-7))
    end

    @testset "Testing propagation correctness" begin
        # Analytical result using complex q factor
        q(q0::Complex, M::Matrix) = (M[1] * q0 + M[3]) / (M[2] * q0 + M[4])
        R(q::Complex) = 1/real(1/q)
        w(q::Complex, λ, n=1) = sqrt(-λ / (π * n * imag(1 / q)))
        propagate_ABCD(d) = [1 d; 0 1]
        lensmaker_ABCD(f) = [1 0; -1/f 1]
        # Beam parameters
        λ = 1000e-9
        w0 = 1e-3
        M2 = 1
        zr = SCDI.rayleigh_range(λ, w0, M2)
        # Lens parameters
        R1 = 1
        R2 = 1
        lens_y_location = 0.1
        nl = 1.5
        f = SCDI.lensmakers_eq(R1, -R2, nl)
        # Stuff
        dy = 0.001
        ys = 0:dy:1.5
        w_analytical = Vector{Float64}(undef, length(ys))
        R_analytical = Vector{Float64}(undef, length(ys))
        # Propagate using ABCD formalism
        q0 = 0 + zr * im
        for i = 1:length(ys)
            w_analytical[i] = w(q0, λ)
            R_analytical[i] = R(q0)
            # catch first lens
            if i * dy == lens_y_location
                q0 = q(q0, lensmaker_ABCD(f))
                continue
            end
            q0 = q(q0, propagate_ABCD(dy))
        end
    
        # Numerical result
        tl = SCDI.ThinLensSDF(R1, R2, 0.025)
        lens = SCDI.Prism(uuid4(), tl, x -> nl)
        system = SCDI.System([lens])
        SCDI.translate3d!(SCDI.shape(lens), [0, lens_y_location, 0])
        # Create and solve beam, calculate beam parameters
        gauss = SCDI.GaussianBeamlet(SCDI.Ray([0, 0, 0], [0, 1, 0]), λ*1e9, w0*1e3, support=[1, 0, 0], M2=1)
        SCDI.solve_system!(system, gauss)
        w_numerical, R_numerical, ψ_numerical = SCDI.gauss_parameters(gauss, ys)
        # Compare beam radius to within 1 μm
        @test all(isapprox.(w_analytical, w_numerical, atol=1e-6))
        # Compare if radius of curvature agreement to within 1 cm above 95%
        temp = isapprox.(R_analytical, R_numerical, atol=1e-2)
        @test sum(temp)/length(temp) > 0.95
        # Compare if Gouy phase zero at waists
        @test isapprox(ψ_numerical[1], 0, atol=1e-3)
        @test isapprox(ψ_numerical[findmin(w_numerical)[2]], 0, atol=1e-3)
    end
end