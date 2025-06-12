using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics
using AbstractTrees

const BMO = BeamletOptics

@testset "Utilities" begin
    @testset "Testing normal3d" begin
        @testset "normal3d(v) with Array" begin
            v = [1.0, 1, 1]
            k = BMO.normal3d(v)
            @test dot(v, k)≈0 atol=1e-14
            @test norm(k) ≈ 1
        end

        @testset "normal3d(v) with Point3" begin
            v = Point3(1.0, 1, 1)
            k = BMO.normal3d(v)
            @test dot(v, k)≈0 atol=1e-14
            @test norm(k) ≈ 1
        end

        @testset "normal3d(v, w) right hand rule and unit length" begin
            orth = BMO.normal3d([2, 0, 0], [0, 0, 1])
            @test isapprox(orth, [0, -1, 0])
        end
    end

    @testset "Testing rotate3d for clockwise dir. and conservation of length" begin
        Rot = BMO.rotate3d([0, 0, 1], π / 2)
        @test isapprox(Rot * [1, 0, 0], [0, 1, 0])
    end

    @testset "Testing align3d for rotation and conservation of length" begin
        # Start vector must have unit length!
        start = [1, 0, 0]
        # Test parallel case
        target = [1, 0, 0]
        T = BMO.align3d(start, target)
        @test T * start ≈ target
        # Test parallel opposite case
        target = [-1, 0, 0]
        T = BMO.align3d(start, target)
        @test T * start ≈ target
        # Test norm and 45° rotation
        target = [1.0, 1.0, 0.0]
        T = BMO.align3d(start, target)
        @test T * start ≈ normalize(target)
    end

    @testset "Testing angle3d for resulting angle" begin
        a = BMO.angle3d([1, 0, 0], [0, 0, 1])
        @test isapprox(a, π / 2)
    end

    @testset "Testing line_point_distance3d and isinfrontof" begin
        pos = [0, 0, 0]
        dir = [1, 0, 0]
        point = [5, 1, 1]
        d = BMO.line_point_distance3d(pos, dir, point)
        @test isapprox(d, √2)

        @test BMO.isinfrontof(point, pos, dir) == true
    end

    @testset "Testing reflection3d" begin
        for dx in -1:1, dy in -1:1
            @test isapprox(
                BMO.reflection3d([dx, dy, 1], [0, 0, -1]), [dx, dy, -1])
        end
    end

    @testset "Testing refraction3d" begin
        normal = [0, 0, 1]
        # Define angle range
        small_angles = 0:1e-7:5e-5
        large_angles = last(small_angles):(π / 1000):(π / 2)
        θs = cat(small_angles, large_angles; dims=1)
        @testset "Test from vacuum into medium" begin
            n1 = 1.0
            n2 = 1.62286
            θ_num = similar(θs)
            θ_ana = similar(θs)
            TIR = Vector{Bool}(undef, length(θs))
            for (i, θ1) in enumerate(θs)
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out, TIR[i] = BMO.refraction3d(dir_in, normal, n1, n2)
                θ_num[i] = BMO.angle3d(-normal, dir_out)
                # 2D-equation for refraction validation
                θ_ana[i] = asin(n1 / n2 * sin(θ1))
        
            end
            @test isapprox(θ_num, θ_ana)
            @test all(TIR .== false)
        end
        @testset "Test from medium into vacuum" begin
            n1 = 1.62286
            n2 = 1.0
            θ2 = similar(θs)
            θ3 = similar(θs)
            TIR_num = Vector{Bool}(undef, length(θs))
            TIR_ana = similar(TIR_num)
            for (i, θ1) in enumerate(θs)
                dir_in = [sin(θ1), 0, -cos(θ1)]
                dir_out, TIR_num[i] = BMO.refraction3d(dir_in, normal, n1, n2)
                if θ1 > asin(n2 / n1)
                    # Test for total reflection
                    θ2[i] = BMO.angle3d(dir_out, normal)
                    θ3[i] = θ1
                    TIR_ana[i] = true
                else
                    # Test for refraction
                    θ2[i] = BMO.angle3d(-normal, dir_out)
                    θ3[i] = asin(n1 / n2 * sin(θ1))
                    TIR_ana[i] = false
                end
            end
            @test isapprox(θ2, θ3)
            @test TIR_ana == TIR_num
        end
    end

    @testset "Fresnel equations" begin
        @testset "Vacuum-glass: normal incidence" begin
            n = 1.5
            θ = 0.0
            rs, rp, ts, tp = BMO.fresnel_coefficients(θ, n)
            @test real(rs) ≈ (1 - n) / (1 + n)
            @test real(rs) ≈ real(rp)
            @test real(tp) ≈ 2 / (1 + n)
            @test real(tp) ≈ real(ts)
        end

        @testset "Vacuum-glass: Brewster angle" begin
            n = 1.5
            θb = atan(n)
            rs, rp, ts, tp = BMO.fresnel_coefficients(θb, n)
            @test real(rp) ≈ 0
        end

        @testset "Vacuum-glass: grazing incidence" begin
            n = 1.5
            θ = π / 2
            rs, rp, ts, tp = BMO.fresnel_coefficients(θ, n)
            @test real(rs) ≈ -1
            @test real(rp) ≈ 1
            @test real(ts) ≈ 0
            @test real(tp)≈0 atol=2e-16
        end

        @testset "Glass-vacuum: normal incidence" begin
            n = 1 / 1.5
            θ = 0.0
            rs, rp, ts, tp = BMO.fresnel_coefficients(θ, n)
            @test real(rs) ≈ (1 - n) / (1 + n)
            @test real(rs) ≈ real(rp)
            @test real(tp) ≈ 2 / (1 + n)
            @test real(tp) ≈ real(ts)
        end

        @testset "Glass-vacuum: Brewster angle" begin
            n = 1 / 1.5
            θb = atan(n)
            rs, rp, ts, tp = BMO.fresnel_coefficients(θb, n)
            @test real(rp)≈0 atol=2e-16
        end

        @testset "Glass-vacuum: Total internal reflection" begin
            n = 1 / 1.5
            θc = asin(n)
            rs, rp, ts, tp = BMO.fresnel_coefficients(θc, n)
            @test BMO.is_internally_reflected(rp, rs)
            @test real(rs) ≈ 1
            @test real(rp) ≈ -1
            @test real(ts) ≈ 2
            @test real(tp)≈3 atol=1e-15
        end
    end

    @testset "Ref. index utils" begin
        # Test data (NLAK22)
        T1 = Float32
        T2 = Float64
        lambdas = T1.([488e-9, 707e-9, 1064e-9])
        indices = T2.([1.6591, 1.6456, 1.6374])
        ref_index = DiscreteRefractiveIndex(lambdas, indices)

        @testset "DiscreteRefractiveIndex" begin
            @test isdefined(BMO, :DiscreteRefractiveIndex)
            @test isa(ref_index, DiscreteRefractiveIndex{T2})
            @test ref_index(lambdas[2]) == indices[2]
            @test_throws KeyError ref_index(lambdas[1] + 1e-9)
            # Test constructor
            @test_throws ArgumentError DiscreteRefractiveIndex([1], [1, 2])
        end

        @testset "Test ref. helper function" begin
            f1(x::Float64) = x              # fail
            f2(x::Union{Int, Float64}) = x  # fail
            f3(x::Real) = "a"               # fail
            f4(x) = x, 1                    # fail
            f5(x) = x                       # pass
            # Test if illegal functions are detected
            @test_throws ArgumentError BMO.test_refractive_index_function(f1)
            @test_throws ArgumentError BMO.test_refractive_index_function(f2)
            @test_throws ArgumentError BMO.test_refractive_index_function(f3)
            @test_throws ArgumentError BMO.test_refractive_index_function(f4)
            @test isnothing(BMO.test_refractive_index_function(f5))
            @test isnothing(BMO.test_refractive_index_function(ref_index))
        end
    end

    @testset "Numerical aperture" begin
        n0 = 1.5
        theta = deg2rad(45)
        BMO.numerical_aperture(theta, n0) == n0 * sin(theta)
    end
end

@testset "Types" begin
    @test isdefined(BMO, :AbstractShape)
    @test isdefined(BMO, :AbstractObject)
    @test isdefined(BMO, :AbstractObjectGroup)
    @test isdefined(BMO, :AbstractRay)
    @test isdefined(BMO, :AbstractBeam)
    @test isdefined(BMO, :AbstractSystem)
    @test isdefined(BMO, :Intersection)
    @test isdefined(BMO, :Hint)
    @test isdefined(BMO, :AbstractInteraction)

    # Generate test structs
    struct TestSystem <: BMO.AbstractSystem end

    @testset "AbstractSystem" begin
        # no tests
        sys = TestSystem()
    end

    mutable struct TestRay{T} <: BMO.AbstractRay{T}
        pos::Vector{T}
        dir::Vector{T}
        λ::T
        n::T
    end

    TestRay(pos::AbstractArray{T}, dir::AbstractArray{T}) where {T} = TestRay(
        pos, dir, T(1000e-9), T(1))

    Base.length(::TestRay) = π

    @testset "AbstractRay" begin
        r = TestRay([0.0, 0, 0], [1.0, 0, 0])
        # Test getters
        @test position(r) == r.pos
        @test BMO.direction(r) == r.dir
        @test BMO.wavelength(r) == r.λ
        @test BMO.refractive_index(r) == r.n
        # Test setters
        n_pos = [1, 1, 1]
        n_dir = [1.0, 1, 0]
        n_lam = 532e-9
        n_rfi = 1.5
        BMO.position!(r, n_pos)
        BMO.direction!(r, n_dir)
        BMO.wavelength!(r, n_lam)
        BMO.refractive_index!(r, n_rfi)
        @test position(r) == n_pos
        @test BMO.direction(r) ≈ n_dir .* (sqrt(2) / 2)
        @test BMO.wavelength(r) == n_lam
        @test BMO.refractive_index(r) == n_rfi
        @test length(r) == π
        # Test ray-plane intersection
        plane_pos = [1, 0, -1]
        plane_nml_1 = [-1, 0, 1]
        plane_nml_2 = [1, 0, 0]
        plane_nml_3 = [0, 0, 1]
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        is_1 = BMO.intersect3d(plane_pos, plane_nml_1, ray)
        is_2 = BMO.intersect3d(plane_pos, plane_nml_2, ray)
        is_3 = BMO.intersect3d(plane_pos, plane_nml_3, ray)
        @test length(is_1) == 2
        @test length(is_2) == 1
        @test isnothing(is_3)
    end

    mutable struct TestBeam{T} <: BMO.AbstractBeam{T, TestRay{T}}
        parent::BMO.Nullable{TestBeam}
        children::Vector{TestBeam}
    end

    TestBeam() = TestBeam{Float64}(nothing, Vector{TestBeam{Float64}}())

    @testset "AbstractBeam" begin
        # Create beam tree
        root = TestBeam()
        cb1 = TestBeam()
        cb2 = TestBeam()
        group = [cb1, cb2]
        cb3 = TestBeam()
        # Add children to root, cb2
        BMO.children!(root, group)
        BMO.children!(cb2, cb3)
        # Test tree structure
        @test treeheight(root) == 2
        @test treebreadth(root) == 2
        @test BMO.children(root) == group
        @test first(BMO.children(cb2)) === cb3
        # Test parent connection
        @test AbstractTrees.parent(root) === nothing
        @test AbstractTrees.parent(cb1) === root
        @test AbstractTrees.parent(cb2) === root
        @test AbstractTrees.parent(cb3) === cb2
        # Replace bottom child
        cbr = TestBeam()
        @test_throws "_modify_beam_head not implemented for $(typeof(cb2))" BMO.children!(
            cb2,
            cbr)
        # Test child removal
        BMO._drop_beams!(cb2)
        @test isempty(BMO.children(cb2))
        # Stuff
        @test_throws "_last_beam_intersection not implemented for $(typeof(cb2))" BMO._last_beam_intersection(cb2)
    end

    mutable struct TestShapeless{T} <: BMO.AbstractShape{T}
        pos::Vector{T}
        dir::Matrix{T}
    end

    TestShapeless() = TestShapeless{Float64}(zeros(3), Matrix{Float64}(I, 3, 3))

    @testset "AbstractShape" begin
        pos = zeros(3)
        dir = Matrix{Float64}(I, 3, 3)
        shape = TestShapeless(pos, dir)
        # Test get/set
        @test BMO.position(shape) == pos
        @test BMO.orientation(shape) == dir
        n_pos = [1, 1, 1]
        n_dir = BMO.rotate3d([0, 0, 1], π / 4)
        BMO.position!(shape, n_pos)
        BMO.orientation!(shape, n_dir)
        @test position(shape) == n_pos
        @test BMO.orientation(shape) == n_dir
        # Test translation
        translate3d!(shape, n_pos)
        @test position(shape) == 2 * n_pos
        reset_translation3d!(shape)
        @test position(shape) == zeros(3)
        # Test rotation for counter-clockwise in right-hand coord. system
        dir = Matrix{Float64}(I, 3, 3)
        BMO.orientation!(shape, dir)
        rotate3d!(shape, [0, 0, 1], deg2rad(45))
        @test all(BMO.orientation(shape)[[1, 2, 5]] .≈ sqrt(2) / 2)
        @test BMO.orientation(shape)[4] ≈ -sqrt(2) / 2
        rotate3d!(shape, [0, 0, 1], deg2rad(135))
        @test BMO.orientation(shape)[1:4:9] == [-1, -1, 1]
        xrotate3d!(shape, π)
        yrotate3d!(shape, π)
        zrotate3d!(shape, π)
        reset_rotation3d!(shape)
        @test BMO.orientation(shape) == dir
        # Test align3d
        target_vec = normalize([1, 1, 1])
        align3d!(shape, target_vec)
        @test BMO.orientation(shape)[:, 2] ≈ target_vec
        reset_rotation3d!(shape)
        # The following test are expected to do nothing but not throw exceptions
        ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
        @test_logs (:warn, "No intersect3d method defined for:") BMO.intersect3d(
            shape,
            ray)

        @testset "Testing AbstractRay - AbstractShape" begin
            shape = TestShapeless([1, 0, 0], Matrix{Int}(I, 3, 3))
            ray = TestRay([0.0, 0, 0], [1.0, 0, 0])
            @test BMO.isinfrontof(shape, ray) == true
            BMO.direction!(ray, -[1, 0, 0])
            @test BMO.isinfrontof(shape, ray) == false
            BMO.direction!(ray, [0, 1, 0])
            @test BMO.isinfrontof(shape, ray) == false
            BMO.direction!(ray, [1.0, 1, 0])
            @test BMO.isinfrontof(shape, ray) == true
            @test norm(BMO.direction(ray)) ≈ 1
        end
    end

    struct TestObject{T, S <: BMO.AbstractShape{T}} <:
           BMO.AbstractObject{T, S}
        shape::S
    end

    TestObject() = TestObject(TestShapeless())

    @testset "AbstractObject" begin
        object = TestObject()
        @test isa(BMO.shape(object), TestShapeless)
        # Test forwarding of kin. API to object shape
        @test position(object) ==
              position(BMO.shape(object))
        @test position(object) ==
              position(BMO.shape(object))
        translate3d!(object, ones(3))
        rotate3d!(object, [0, 0, 1], π)
        @test position(object) == ones(3)
        @test BMO.orientation(object)[1:4:9] == [-1, -1, 1]
        reset_translation3d!(object)
        reset_rotation3d!(object)
        @test position(object) == zeros(3)
        @test BMO.orientation(object)[1:4:9] == ones(3)
        # Test translate_to3d
        target_pos = [1, 3, 9]
        translate_to3d!(object, target_pos)
        @test position(object) == target_pos

        @testset "Testing interact3d" begin
            sys = TestSystem()
            obj = TestObject()
            ray = TestRay(zeros(3), ones(3))
            beam = TestBeam()
            @test_logs (:warn, "No interact3d method defined for:") BMO.interact3d(
                sys,
                obj,
                beam,
                ray)===nothing
        end
    end
end

@testset "Rays" begin
    # Testing constructor
    pos = [0, 0, 0]
    dir = [1.0, 1, 0]
    ray = Ray(pos, dir)
    @test ismutable(ray)
    @test isa(ray, Ray{Float64})
    @test isnothing(ray.intersection)
    @test isinf(length(ray))
    @test isapprox(norm(ray.dir), 1)
    # Test helper functions
    @test BMO.line_point_distance3d(ray, [1, 1, 0]) == 0
    @test BMO.line_point_distance3d(ray, [-1, 1, 0]) == sqrt(2)

    @testset "Testing isentering" begin
        r1 = Ray([0, 0, 0], [0, 1, 0])
        r2 = Ray([0, 0, 0], [0, 1, 0])
        r3 = Ray([0, 0, 0], [0, 1, 0])
        i1 = BMO.Intersection(nothing, nothing, 0.0, Point3(0, 1.0, 0))
        i2 = BMO.Intersection(nothing, nothing, 0.0, Point3(0, -1.0, 0))
        BMO.intersection!(r1, i1)
        BMO.intersection!(r2, i2)
        @test !BMO.isentering(r1)
        @test BMO.isentering(r2)
        @test !BMO.isentering(r3)
    end

    @testset "Testing refraction3d" begin
        n1 = 1
        n2 = 1.5
        dir = normalize([0, 1, 0])
        ray = Ray(zeros(3), dir)
        nml = normalize(Point3{Float64}(0, -1, 1))
        BMO.intersection!(
            ray, BMO.Intersection(nothing, nothing, 1.0, nml))
        @test BMO.refraction3d(dir, nml, n1, n2) ==
              BMO.refraction3d(ray, n2)
        # test for correct exit normal flip
        nml *= -1
        BMO.intersection!(
            ray, BMO.Intersection(nothing, nothing, 1.0, nml))
        @test BMO.refraction3d(dir, -nml, n1, n2) ==
              BMO.refraction3d(ray, n2)
    end
end

@testset "Beams" begin
    is = BMO.Intersection(1.0, zeros(3))
    r1 = Ray([0.0, 0, 0], [1, 0, 0])
    r2 = Ray([1.0, 0, 0], [0, 1, 0])
    r3 = Ray([1.0, 1, 0], [0, 0, 1])
    r4 = Ray([1.0, 1, 1], [1, 0, 0])
    BMO.intersection!(r1, is)
    BMO.intersection!(r2, is)
    BMO.intersection!(r3, is)
    BMO.intersection!(r4, nothing)
    # Test beam
    beam = Beam(r1)
    push!(beam, r2)
    push!(beam, r3)
    push!(beam, r4)
    @test length(beam) == 3
    @test BMO.point_on_beam(beam, 0) == ([0, 0, 0], 1)
    @test BMO.point_on_beam(beam, 1) == ([1, 0, 0], 2)
    @test BMO.point_on_beam(beam, 2) == ([1, 1, 0], 3)
    @test BMO.point_on_beam(beam, 3) == ([1, 1, 1], 4)
    @test BMO.point_on_beam(beam, 10) == ([8, 1, 1], 4)
    @test BMO.isparentbeam(beam, r2) == true
end

@testset "Beam groups" begin
    @testset "AbstractBeamGroup definitions" begin
        @test isdefined(BMO, :AbstractBeamGroup)
        
        struct BeamTestGroup{T} <: BMO.AbstractBeamGroup{T, Ray{T}}
            central_beam::BMO.Beam{T, Ray{T}}
        end
        
        struct BeamTestSystem <: BMO.AbstractSystem end
        BMO.objects(::BeamTestSystem) = [nothing]
        BMO.intersect3d(::Nothing, ::BMO.AbstractRay) = nothing
        BMO.interact3d(::Nothing) = nothing
        
        BMO.beams(tg::BeamTestGroup) = [tg.central_beam]
        
        pos = [0,0,0]
        dir = [0,1,0]
        lambda = 1064e-9
        
        tg = BeamTestGroup(BMO.Beam(pos, dir, lambda))
        
        @test position(tg) == pos
        @test BMO.direction(tg) == dir
        @test BMO.wavelength(tg) == lambda
        
        ts = BeamTestSystem()
        
        @test isnothing(solve_system!(ts, tg))
    end

    is_strictly_sorted(v) = all(diff(v) .> 0)

    @testset "Point source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 1]
        alpha = deg2rad(2)
        NA = BMO.numerical_aperture(alpha)
        num_rays = 1000
        num_rings = 10
        
        source = PointSource(pos, dir, alpha, lambda; num_rays, num_rings)
        
        @testset "Testing point source getters" begin
            @test BMO.numerical_aperture(source) == NA
            @test position(source) == pos
            @test BMO.direction(source) ≈ normalize(dir)
        end
        
        @testset "Testing point source max. spread" begin
            last_ray = first(BMO.rays(last(BMO.beams(source))))
            @test BMO.angle3d(dir, BMO.direction(last_ray)) ≈ alpha atol = 1e-14
            @test position(last_ray) == pos
            @test length(BMO.beams(source)) == num_rays
        end
        
        @testset "Testing generated angles" begin
            # Test for center ray 0 deg, generated spread angles etc.
            directions = BMO.direction.(first.(BMO.rays.(BMO.beams(source))))
            angles = BMO.angle3d.(Ref(dir), directions)
            generated_angles = unique(round.(angles, digits=11))
            required_angles = LinRange(0, alpha, num_rings)
            @test generated_angles ≈ required_angles
            # Test for only one center ray, strictly increasing rays per ring
            rays_per_angle = zeros(Int, length(required_angles))
            for (i, _required) in enumerate(required_angles)
                rays_per_angle[i] = count(_required .≈ angles)
            end
            @test first(rays_per_angle) == 1
            @test is_strictly_sorted(rays_per_angle)
            # Test if all generated rays are unique
            @test length(unique(directions)) == length(directions)
        end
    
        @testset "Testing throw errors" begin
            @test_throws ErrorException PointSource(pos, dir, 1.1*π, lambda; num_rays, num_rings)
            @test_throws ErrorException PointSource(pos, dir, alpha, lambda; num_rays=100, num_rings=10)
        end
    end

    @testset "Collimated source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 0]
        diameter = 2BMO.inch
        num_rays = 500
        num_rings = 5
        
        source = CollimatedSource(pos, dir, diameter; num_rays, num_rings)
        
        @testset "Testing coll. source getters" begin 
            @test BMO.diameter(source) == diameter
            @test position(source) == pos
            @test BMO.direction(source) == dir
            @test BMO.wavelength(source) == 1e-6
        end
        
        @testset "Testing coll. source max. spread diameter" begin
            last_ray = first(BMO.rays(last(BMO.beams(source))))
            @test BMO.direction(last_ray) == dir
            @test norm(position(last_ray) - pos) ≈ diameter/2
            @test length(BMO.beams(source)) == num_rays
        end
        
        @testset "Testing coll. source generated positions" begin
            # Test for center ray 0 offset, generated radii
            positions = position.(first.(BMO.rays.(BMO.beams(source))))
            radii = norm.(positions .- Ref(pos))
            generated_pos = unique(round.(radii, digits=11))
            required_pos = LinRange(0, diameter/2, num_rings)
            @test generated_pos ≈ required_pos
            # Test for only one center ray, strictly increasing rays per ring
            rays_per_radius = zeros(Int, length(required_pos))
            for (i, _required) in enumerate(required_pos)
                rays_per_radius[i] = count(_required .≈ radii)
            end
            @test first(rays_per_radius) == 1
            @test is_strictly_sorted(rays_per_radius)
            # Test if all generated rays are unique
            @test length(unique(positions)) == length(positions)
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException CollimatedSource(pos, dir, diameter; num_rays=100, num_rings=10)
        end
    end
end

@testset "Mesh" begin
    # NOTE: the "Mesh" testset is mutating. Errors/fails might lead to subsequent tests failing too!
    @test isdefined(BMO, :AbstractMesh)
    @test isdefined(BMO, :Mesh)

    # Generate cube since types are defined
    foo = BMO.CubeMesh(1) # test cube
    bar = BMO.CubeMesh(1) # reference cube

    to_origin = -0.5 * [1, 1, 1]

    @testset "Testing AbstractMesh getters" begin
        @test typeof(foo) == BMO.Mesh{Float64}
        @test BMO.vertices(foo) == foo.vertices
        @test BMO.faces(foo) == foo.faces
        @test BMO.orientation(foo) == foo.dir
        @test position(foo) == foo.pos
        @test BMO.scale(foo) == foo.scale
    end

    @testset "Testing translate3d!" begin
        translate3d!(foo, to_origin) # move COG to origin
        @test minimum(BMO.vertices(foo)[:, 1]) == -0.5
        @test minimum(BMO.vertices(foo)[:, 2]) == -0.5
        @test minimum(BMO.vertices(foo)[:, 3]) == -0.5
        @test maximum(BMO.vertices(foo)[:, 1]) == 0.5
        @test maximum(BMO.vertices(foo)[:, 2]) == 0.5
        @test maximum(BMO.vertices(foo)[:, 3]) == 0.5
        @test all(position(foo) .== -0.5)
    end

    @testset "Testing set_new_origin3d!" begin
        BMO.set_new_origin3d!(foo)
        @test position(foo) == zeros(3)
    end

    @testset "Testing x/y/zrotate3d!" begin
        @testset "Testing rotate3d!" begin
            rotate3d!(foo, [1, 0, 0], π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), √2 / 2)
            # Return to original rotation
            rotate3d!(foo, [1, 0, 0], -π / 4)
        end

        @testset "Testing xrotate3d!" begin
            xrotate3d!(foo, π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -√2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), √2 / 2)
        end

        @testset "Testing yrotate3d!" begin
            yrotate3d!(foo, π / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -√2 / 2)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), √2 / 2)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 0.5)
        end

        @testset "Testing zrotate3d!" begin
            zrotate3d!(foo, π / 4)
            @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -0.5)
            @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 2]), 0.5)
            @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 0.5)
        end

        # Testing orientation of dir matrix
        @test BMO.orientation(foo)[[3, 5, 7]] == [-1, 1, 1]
    end

    # center bar reference cube at origin
    translate3d!(bar, to_origin)
    BMO.set_new_origin3d!(bar)

    @testset "Testing reset_rotation3d!" begin
        translate3d!(foo, [1, 2, 3])
        reset_translation3d!(foo)
        reset_rotation3d!(foo)
        @test position(foo) == zeros(3)
        @test BMO.orientation(foo) ≈ BMO.orientation(bar)
        @test BMO.vertices(foo) ≈ BMO.vertices(bar)
    end

    @testset "Testing align3d!" begin
        align3d!(foo, normalize([0, 1, 1]))
        @test position(foo) == zeros(3)
        @test BMO.orientation(foo)[:, 1] ≈ [1, 0, 0]
        @test BMO.orientation(foo)[:, 2] ≈ [0, √2 / 2, √2 / 2]
        @test BMO.orientation(foo)[:, 3] ≈ [0, -√2 / 2, √2 / 2]
        reset_rotation3d!(foo)
    end

    @testset "Testing normal" begin
        normal = BMO.normal3d(foo, 1)
        @test isapprox(normal, [0, 0, -1])
    end

    @testset "Testing scale3d!" begin
        BMO.scale3d!(foo, 2)
        @test isapprox(minimum(BMO.vertices(foo)[:, 1]), -1)
        @test isapprox(minimum(BMO.vertices(foo)[:, 2]), -1)
        @test isapprox(minimum(BMO.vertices(foo)[:, 3]), -1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 1]), 1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 2]), 1)
        @test isapprox(maximum(BMO.vertices(foo)[:, 3]), 1)
        @test BMO.scale(foo) == 2
    end

    @testset "Testing Moeller-Trumbore algorithm" begin
        t = 5
        face = [1 1 t
                -1 1 t
                0 -1 t]
        # ray at origin pointing along z-axis
        pos = [0.0, 0, 0]
        dir = [0.0, 0, 1]
        ray = Ray(pos, dir)
        # Preallocate memory
        @test isapprox(BMO.MoellerTrumboreAlgorithm(face, ray), t)
        # Check allocations (WARNING: function must have been compiled once for before this test!)
        alloc = @allocated BMO.MoellerTrumboreAlgorithm(face, ray)
        if alloc > 16
            @warn "Allocated number of bytes for MTA larger than expected!" alloc
        end
    end
    @testset "Testing intersect3d" begin
        # Setup test cube and ray
        cube = BMO.CubeMesh(1)
        translate3d!(cube, -0.5 * [1, 1, 1])
        BMO.set_new_origin3d!(cube)
        ray_pos = zeros(3)
        ray_dir = [1.0, 0, 0]
        ray = Ray(ray_pos, ray_dir)
        # Rotate cube 360°, calculate intersection distance
        θ = 0:1:359
        l = zeros(length(θ))
        for (i, ~) in enumerate(θ)
            intersection = BMO.intersect3d(cube, ray)
            l[i] = length(intersection)
            zrotate3d!(cube, deg2rad(step(θ)))
        end
        # Test if 0/45° distances are correct
        @test all(l[1:90:end] .≈ BMO.scale(cube) * 1 / 2)
        @test all(l[(1 + 45):90:end] .≈ BMO.scale(cube) * sqrt(2) / 2)
    end

    @testset "Testing intersect3d - part 2" begin
        t = 5
        s = 1 # scale/2
        cube = BMO.CubeMesh(2 * s)
        # Move cube COG to origin
        translate3d!(cube, -[s, s, s])
        BMO.set_new_origin3d!(cube)
        # Align cube edge at t units from origin
        translate3d!(cube, [t + s, 0, 0])
        pos = [0, 0, 0]
        steps = 10
        for z in (-s):(s / steps):s
            # Ray constructed each time for unit-length dir
            dir = [t, 0, z]
            ray = Ray(pos, dir)
            @test isapprox(BMO.intersect3d(cube, ray).t, sqrt(t^2 + z^2))
        end
    end

    @testset "Testing constructors" begin
        @testset "Testing RectangularFlatMesh" begin
            rfm = BMO.RectangularFlatMesh(2.0, 1)
            @test BMO.vertices(rfm) == [1 0 0.5; 1 0 -0.5; -1 0 -0.5; -1 0 0.5]
            @test BMO.normal3d(rfm, 1) == [0, 1, 0]
        end

        @testset "Testing QuadraticFlatMesh" begin
            qfm = BMO.QuadraticFlatMesh(4.0)
            @test BMO.vertices(qfm) == [2 0 2; 2 0 -2; -2 0 -2; -2 0 2]
            @test BMO.normal3d(qfm, 1) == [0, 1, 0]
        end
    end

    @testset "Testing CircularFlatMesh" begin
        # Testing constructor
        n = 4
        cm = BMO.CircularFlatMesh(1.0f0, n)
        v = BMO.vertices(cm)
        f = BMO.faces(cm)

        # testing vertices
        @test v[1, :] ≈ zeros(3)
        @test v[2, :] ≈ [1, 0, 0]
        @test v[3, :] ≈ [0, 0, 1]
        @test v[4, :] ≈ [-1, 0, 0]
        @test v[5, :] ≈ [0, 0, -1]

        # testing faces
        @test f[:, 1] == ones(4)
        @test f[:, 2] == [2, 3, 4, 5]
        @test f[:, 3] == [3, 4, 5, 2]

        # testing normal vectors
        for i in 1:n
            @test BMO.normal3d(cm, i) ≈ [0, -1, 0]
        end
    end
end

@testset "SDFs" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

    # Orientation-less test point sdf
    mutable struct TestPointSDF{T} <: BMO.AbstractSDF{T}
        position::Point3{T}
        orientation::Matrix{T}
    end

    TestPointSDF(p::AbstractArray{T}) where {T} = TestPointSDF{T}(
        Point3{T}(p), Matrix{T}(I, 3, 3))
    TestPointSDF(T = Float64) = TestPointSDF{T}(Point3{T}(0), Matrix{T}(I, 3, 3))

    BMO.position(tps::TestPointSDF) = tps.position
    BMO.position!(tps::TestPointSDF{T}, new::Point3{T}) where {T} = (tps.position = new)

    BMO.orientation(tps::TestPointSDF) = tps.orientation
    BMO.orientation!(tps::TestPointSDF{T}, new::Matrix{T}) where {T} = (tps.orientation = new)

    BMO.transposed_orientation(tps::TestPointSDF) = transpose(tps.orientation)
    BMO.transposed_orientation!(::TestPointSDF, ::Any) = nothing

    orientation(::TestPointSDF{T}) where {T} = Matrix{T}(I, 3, 3)
    orientation!(::TestPointSDF, ::Any) = nothing

    function BMO.sdf(tps::TestPointSDF, point)
        p = BMO._world_to_sdf(tps, point)
        return norm(p)
    end

    @testset "Testing kinematics and transforms" begin
        point = TestPointSDF()
        t = 10
        θ = deg2rad(30)
        translate3d!(point, [t, 0, 0])
        rotate3d!(point, [0, 1, 0], θ)
        pt = BMO._world_to_sdf(point, [0, 0, 0])
        @test pt[1] ≈ -t * cos(θ)
        @test pt[2] ≈ 0
        @test pt[3] ≈ -t * sin(θ)
    end

    @testset "Testing intersect3d" begin
        t = 10.0
        point = TestPointSDF(zeros(3))
        translate3d!(point, [t, 0, 0])

        r1 = Ray(zeros(3), [1.0, 0, 0])
        r2 = Ray(zeros(3), [1.0, 1, 0])
        r3 = Ray(zeros(3), [1.0, 0, 1])

        i1 = BMO.intersect3d(point, r1)
        i2 = BMO.intersect3d(point, r2)
        i3 = BMO.intersect3d(point, r3)

        @test length(i1) == t
        @test isnothing(i2)
        @test isnothing(i3)
    end

    @testset "Testing normal3d" begin
        point = TestPointSDF(zeros(3))
        offset = [5, 0, 0]
        translate3d!(point, offset)
        p1 = [1, 0, 0]
        p2 = [0, 1, 0]
        p3 = [0, 0, 1]
        @test BMO.normal3d(point, p1 + offset) == p1
        @test BMO.normal3d(point, p2 + offset) == p2
        @test BMO.normal3d(point, p3 + offset) == p3
    end
end

@testset "System" begin
    @testset "Testing implementation" begin
        struct SystemTestBeam{T} <: BMO.AbstractBeam{T, Ray{T}} end
        struct SystemTestObject{T, S} <: BMO.AbstractObject{T, S} end
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
end

@testset "Object groups" begin
    mutable struct TestPoint{T} <: BMO.AbstractShape{T}
        pos::Point3{T}
        dir::Matrix{T}
    end

    TestPoint(position::AbstractArray{T}) where {T <: Real} = TestPoint{T}(
        Point3{T}(position),
        Matrix{T}(I, 3, 3))

    struct GroupTestObject{T <: Real, S <: BMO.AbstractShape{T}} <:
           BMO.AbstractObject{T, S}
        shape::S
    end

    GroupTestObject(position::AbstractArray) = GroupTestObject(TestPoint(position))

    n = 8
    xs = [cos(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]
    ys = [sin(x) for x in LinRange(0, 2pi * (n - 1) / n, n)]

    # Test center with Float32, rest with Float64
    center = GroupTestObject(zeros(Float32, 3))
    circle = ObjectGroup([GroupTestObject([xs[i], ys[i], 0]) for i in eachindex(xs)])

    objects = ObjectGroup([center, circle])

    # Translation test
    target = [3, 0, 0]
    translate_to3d!(objects, target)

    @testset "translate3d" begin
        # Test if all objects/subgroups have been translated
        @test position(objects) == target
        @test position(center) == target
        @test position(circle) == target
        for (i, obj) in enumerate(BMO.objects(circle))
            @test position(obj) == [xs[i], ys[i], 0] + target
        end
    end

    # Rotation test
    angle = 2π / n
    rotate3d!(objects, [0, 0, 1], angle)

    @testset "rotate3d" begin
        # Test if all objects/subgroups have been rotated relative to the origin
        Rt = BMO.rotate3d([0, 0, 1], angle)
        xt = circshift(xs, -1)
        yt = circshift(ys, -1)
        @test BMO.orientation(objects) == Rt
        @test BMO.orientation(center) ≈ Rt
        @test BMO.orientation(circle) == Rt
        for (i, obj) in enumerate(BMO.objects(circle))
            @test BMO.orientation(obj) == Rt
            @test position(obj) ≈ [xt[i], yt[i], 0] + target
        end
    end

    # Reset test
    reset_translation3d!(objects)
    reset_rotation3d!(objects)

    @testset "reset functions" begin
        Ri = Matrix{Float64}(I, 3, 3)
        # Test if objects are reset correctly to initial positioning
        @test position(objects) == zeros(3)
        @test position(center) == zeros(3)
        @test position(circle) == zeros(3)
        @test BMO.orientation(objects) == Ri
        @test BMO.orientation(center) ≈ Ri
        @test BMO.orientation(circle) ≈ Ri
        for (i, obj) in enumerate(Leaves(BMO.objects(circle)))
            @test isapprox(position(obj)[1], xs[i], atol = 5e-16)
            @test isapprox(position(obj)[2], ys[i], atol = 5e-16)
        end
    end

    @testset "System compatibility" begin
        # Test if objects in ObjectGroup are exposed correctly when iterating
        system = System(objects)
        ctr = 0
        # Only the objects within the groups should be exposed
        for obj in BMO.objects(system)
            @test isa(obj, GroupTestObject)
            ctr += 1
        end
        @test ctr == n + 1
    end
end

@testset "Spherical Lenses" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :AbstractSDF)
        @test isdefined(BMO, :SphereSDF)
        @test isdefined(BMO, :CylinderSDF)
        @test isdefined(BMO, :CutSphereSDF)
        @test isdefined(BMO, :ThinLensSDF)
    end

    @testset "Thin lens focal length" begin
        # define thin lens
        R1 = 1
        R2 = 1
        nl = 1.5
        tl = BMO.ThinLensSDF(R1, R2, 0.1)
        translate3d!(tl, [0, -thickness(tl) / 2, 0])
        p = Lens(tl, x -> 1.5)
        system = System(p)

        # compare numerical and analytical focal length
        f_analytical = BMO.lensmakers_eq(R1, -R2, nl)
        zs = -0.04:0.01:0.04
        for (i, z) in enumerate(zs)
            # skip optical axis ray
            if z ≈ 0
                continue
            end
            xs = 0.1:0.1:1.5
            df = zeros(Float64, length(xs))
            ray = Ray([0, -0.5, z], [0, 1, 0], 1e3)
            beam = Beam(ray)
            solve_system!(system, beam)
            # test if numerical and analytical focal length agree
            for (i, x) in enumerate(xs)
                df[i] = BMO.line_point_distance3d(beam.rays[end], [0, x, 0])
            end
            @test xs[findmin(df)[2]] ≈ f_analytical
        end
    end

    @testset "Testing lens constructor" begin
        # Test against Thorlab spherical lenses
        r1 = 34.9e-3
        r2 = -r1
        l = 6.8e-3
        LB1811 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LB1811)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LB1811)) == l
        r1 = Inf
        r2 = -15.5e-3
        l = 8.6e-3
        LA1805 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LA1805)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LA1805)) == l
        r1 = -52.0e-3
        r2 = -r1
        l = 3e-3
        LD1464 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LD1464)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LD1464)) == l
        r1 = Inf
        r2 = 25.7e-3
        l = 3.5e-3
        LC1715 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LC1715)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LC1715)) == l
        r1 = -82.2e-3
        r2 = -32.1e-3
        l = 3.6e-3
        LE1234 = SphericalLens(r1, r2, l)
        @test typeof(BMO.shape(LE1234)) <: BMO.UnionSDF
        @test thickness(BMO.shape(LE1234)) == l
    end

    """Test coma for rotated and translated optical system"""
    function test_coma(ray::BMO.AbstractRay, f0::AbstractArray,
            dir::AbstractArray; atol = 7e-5)
        is = BMO.intersect3d(f0, dir, ray)
        p0 = position(ray) + length(is) * BMO.direction(ray)
        dz = norm(p0 - f0)
        if dz ≤ atol
            return true
        else
            return error("Coma dz=$dz larger than atol=$atol")
        end
    end

    @testset "Testing doublet lenses" begin
        # Define refractive index functions
        λs = [488e-9, 707e-9, 1064e-9]
        NLAK22 = DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
        NSF10 = DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

        function test_doublet(λ, bfl, δf)
            # Thorlabs lens from https://www.thorlabs.com/thorproduct.cfm?partnumber=AC254-150-AB
            AC254_150_AB = SphericalDoubletLens(
                87.9e-3, -105.6e-3, Inf, 6e-3, 3e-3, BMO.inch, NLAK22, NSF10)
            # Rotate and translate to test lens kinematics
            translate3d!(AC254_150_AB, [0.05, 0.05, 0.05])
            xrotate3d!(AC254_150_AB, deg2rad(-60))
            zrotate3d!(AC254_150_AB, deg2rad(45))
            # Define system
            system = System([AC254_150_AB])
            # Define semi-diameter for lens ray bundle, selected for min. spherical aberrations
            z0 = 5e-3
            zs = LinRange(-z0, z0, 30)
            fs = similar(zs)
            # Beam spawn point
            dir = -BMO.orientation(AC254_150_AB.back.shape)[:, 2]       # rotated collimated ray direction
            pos = position(AC254_150_AB.front.shape) + 0.05 * dir  # rotated collimated ray position
            nv = BMO.normal3d(dir)                                     # orthogonal to moved system optical axis
            beam = Beam(pos, -dir, λ)
            # Calculate equivalent back focal length point
            f_z = thickness(AC254_150_AB) + bfl + δf
            f0 = position(AC254_150_AB.front.shape) + f_z * -dir
            for (i, z) in enumerate(zs)
                beam.rays[1].pos = pos + z * nv
                solve_system!(system, beam)
                @test length(BMO.rays(beam)) == 4
                @test BMO.refractive_index.(beam.rays) ==
                      [1, NLAK22(λ), NSF10(λ), 1]
                fs[i] = test_coma(last(BMO.rays(beam)), f0, dir, atol = 1e-6)
            end
            # Test center ray normal vectors
            beam.rays[1].pos = pos + 0 * nv
            solve_system!(system, beam)
            for i in 1:(length(beam.rays) - 1)
                @test abs(dot(beam.rays[i].intersection.n, beam.rays[i].dir)) ≈ 1
            end
            return true
        end
        # Run tests for AC254_150_AB against plot data at https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=12767
        @test test_doublet(488e-9, 143.68e-3, -2.064e-4)
        @test test_doublet(707e-9, 143.68e-3, 0)
        @test test_doublet(1064e-9, 143.68e-3, +7.466e-4)
    end
end

@testset "Surfaces" begin
    function working_distance(lens, offset_z)
        ray = Ray([0.0, -1.0, offset_z], [0.0, 1.0, 0])
        system = System([lens])
        beam = Beam(ray)
        solve_system!(system, beam)

        dist = -beam.rays[end].pos[3] / beam.rays[end].dir[3]
        α = asind(beam.rays[end].dir[3])
        wd = cosd(α) * dist

        return wd
    end

    @testset "Lens construction from surfaces" begin
        ## Thorlabs LA4052, plano-convex
        r1 = 16.1e-3
        r2 = Inf
        l = 8.2e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            l,
            n -> 1.458
        )

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test thickness(BMO.shape(lens).sdfs[1])≈2e-3 atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈29.5e-3 atol=1e-4

        ## Thorlabs LB1761, bi-convex
        r1 = 24.5e-3
        r2 = -24.5e-3
        l = 9.0e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test thickness(shape.sdfs[1])≈1.9e-3 atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈22.2e-3 atol=1e-3

        ## Thorlabs LC1715, plano-concave
        r1 = Inf
        r2 = 25.7e-3
        l = 3.5e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l
        # test edge thickness
        @test BMO.sag(shape.sdfs[2]) +
              thickness(shape.sdfs[1])≈0.006858 atol=1e-4

        ## Thorlabs LD2297, bi-concave
        r1 = -39.6e-3
        r2 = 39.6e-3
        l = 3.0e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        @test BMO.sag(shape.sdfs[2]) + BMO.sag(shape.sdfs[3]) +
              thickness(shape.sdfs[1])≈0.0072 atol=1e-4

        ## Thorlabs LBF254-040, best-form
        r1 = 134.6e-3
        r2 = -24.0e-3
        l = 6.5e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1])≈2.286e-3 atol=1e-4

        ## Thorlabs LE1234, positive meniscus
        r1 = -82.2e-3
        r2 = -32.1e-3
        l = 3.6e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1]) + BMO.sag(shape.sdfs[2])≈2e-3 atol=1e-4

        ## Thorlabs LF1822, negative meniscus
        r1 = -33.7e-3
        r2 = -100.0e-3
        l = 3.0e-3
        d = 25.4e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test edge thickness
        @test thickness(shape.sdfs[1]) +
              BMO.sag(shape.sdfs[2])≈4.7e-3 atol=1e-4

        ## Generic "true" meniscus
        r1 = 103.4371e-3
        r2 = 61.14925e-3
        l = 1.5e-3
        d = 55e-3
        lens = Lens(
            SphericalSurface(r1, d),
            SphericalSurface(r2, d),
            l,
            n -> 1.517
        )
        shape = lens.shape

        # test lens thickness
        @test thickness(lens) ≈ l

        # test ring generation
        NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])
        mm = 1e-3

        s1 = Lens(
            SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm),
            SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm),
            0.5mm,
            NBK7
        )

        @test 2*s1.shape.sdfs[5].hthickness ≈ 0.001134 atol=1e-6

        s2 = Lens(
            SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm),
            SphericalSurface(-5.020mm, 2*2.380mm, 2*2.380mm),
            2.5mm,
            NBK7
        )

        @test 2*s2.shape.sdfs[4].hthickness ≈ 0.001221590 atol=1e-8

        # doublet test case
        s1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
        s2 = SphericalSurface(-3.642mm, 2*3mm)
        s3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)

        dl21 = Lens(s1, s2, 3.4mm, NBK7)
        dl22 = Lens(s2, s3, 1.0mm, NBK7)

        @test 2*dl21.shape.sdfs[4].hthickness ≈ 0.001294398 atol=1e-6
        @test 2*dl22.shape.sdfs[4].hthickness ≈ 0.000723025 atol=1e-6
    end
end

@testset "Aspherical Lenses" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :ConvexAsphericalSurfaceSDF)
        @test isdefined(BMO, :ConcaveAsphericalSurfaceSDF)
    end

    @testset "Testing aspherical lens SDFs" begin
        # Test that the SDF in combination with ray-marching correctly approximates the
        # aspheric surface. Let's use the Thorlabs AL50100J lens for this test case.

        # radius
        R = 50.3583e-3
        # conic constant
        k = -0.789119
        # even aspheric coefficients
        A = [0, 2.10405e-7 * (1e3)^3, 1.76468e-11 * (1e3)^5, 1.02641e-15 * (1e3)^7]
        # center thickness
        ct = 10.2e-3
        # diameter
        d = 50e-3
        # refractive index of BK-7 @ 1310 nm (design wavelength)
        n = 1.5036

        lens = Lens(
            EvenAsphericalSurface(R, d, k, A),
            ct,
            x -> n
        )

        system = System(lens)

        surf_errors = zeros(100)

        for (i, z) in enumerate(range(-0.02, 0.02, 100))
            ray = Ray(Point3(0.0, -0.1, z), Point3(0.0, 1.0, 0))
            beam = Beam(ray)
            solve_system!(system, beam, r_max = 40)

            surf_errors[i] = (position(beam.rays[begin]) + length(beam.rays[begin]) .* BMO.direction(beam.rays[begin]))[2] -
                             BMO.aspheric_equation(ray.pos[3], 1 / R, k, A)
        end

        # FIXME: The atol is actually derived from the raymarching epsilon. If this is puts
        # into a configurable option, this should be changed as well.
        @test all(x -> isapprox(x, 0.0; atol = 1e-10), surf_errors)

        # test if the working distance is correct
        ray = Ray([0.0, -0.1, 0.02], [0.0, 1.0, 0])
        beam = Beam(ray)
        solve_system!(system, beam, r_max = 40)

        dist = -beam.rays[end].pos[3] / beam.rays[end].dir[3]
        α = asind(beam.rays[end].dir[3])
        wd = cosd(α) * dist

        @test wd≈93.2e-3 atol=1e-4
    end

    @testset "Complex aspherical imaging system" begin
        # setup system
        L1 = Lens(
            EvenAsphericalSurface(
                1.054e-3, # r
                1.333024e-3, # d
                -0.14294, # conic
                [0, 0.038162 * (1e3)^3, 0.06317 * (1e3)^5,
                    -0.020792 * (1e3)^7, 0.18432 * (1e3)^9,
                    -0.04827 * (1e3)^11, 0.094529 * (1e3)^13] # coeffs
            ),
            EvenAsphericalSurface(
                2.027e-3, # r
                1.216472e-3, # d
                8.0226, # conic
                [0, 0.0074974 * (1e3)^3, 0.064686 * (1e3)^5,
                    0.19354 * (1e3)^7, -0.50703 * (1e3)^9,
                    -0.34529 * (1e3)^11, 5.9938 * (1e3)^13] # coeffs
            ),
            0.72e-3, # center thickness
            n -> 1.580200
        )

        L2 = Lens(
            EvenAsphericalSurface(
                -3.116e-3, # r
                1.4e-3, # d
                -49.984, # conic
                [0, -0.31608 * (1e3)^3, 0.34755 * (1e3)^5,
                    -0.17102 * (1e3)^7, -0.41506 * (1e3)^9,
                    -1.342 * (1e3)^11, 5.0594 * (1e3)^13, -2.7483 * (1e3)^15] # coeffs
            ),
            EvenAsphericalSurface(
                -4.835e-3, # r
                1.9e-3, # d
                1.6674, # conic
                [0, -0.079727 * (1e3)^3, 0.13899 * (1e3)^5, -0.044057 * (1e3)^7,
                    -0.019369 * (1e3)^9, 0.016993 * (1e3)^11, 0.093716 * (1e3)^13,
                    -0.080329 * (1e3)^15] # coeffs
            ),
            0.55e-3, # center_thickness
            n -> 1.804700
        )

        translate3d!(L2, [0, thickness(L1) + 0.39e-3, 0])

        L3 = Lens(
            EvenAsphericalSurface(
                3.618e-3, # r
                3.04e-3, # d
                -44.874, # conic
                [0, -0.14756 * (1e3)^3, 0.035194 * (1e3)^5, -0.0032262 * (1e3)^7,
                    0.0018592 * (1e3)^9, 0.00036658 * (1e3)^11, -0.00016039 * (1e3)^13,
                    -3.1846e-5 * (1e3)^15] # coeffs
            ),
            EvenAsphericalSurface(
                2.161e-3, # r
                3.7e-3, # d
                -10.719, # conic
                [0, -0.096568 * (1e3)^3, 0.026771 * (1e3)^5, -0.011261 * (1e3)^7,
                    0.0019879 * (1e3)^9, 0.00015579 * (1e3)^11, -0.00012433 * (1e3)^13,
                    1.5264e-5 * (1e3)^15] # coeffs
            ),
            0.7e-3, # center_thickness
            n -> 1.580200
        )

        translate_to3d!(L3, position(L2))
        translate3d!(L3, [0, thickness(L2) + 0.63e-3, 0])

        Filt = Lens(
            CircularFlatSurface(4.2e-3),
            0.15e-3,
            n -> 1.516800
        )

        translate_to3d!(Filt, position(L3))
        translate3d!(Filt, [0, thickness(L3) + 0.19e-3, 0])

        Cover = Lens(
            CircularFlatSurface(4.9e-3),
            0.5e-3,
            n -> 1.469200
        )
        translate_to3d!(Cover, position(Filt))
        translate3d!(Cover, [0, thickness(Filt) + 0.18e-3, 0])

        # test thickness
        @test thickness(L1) ≈ 0.72e-3
        @test thickness(L2) ≈ 0.55e-3
        @test thickness(L3) ≈ 0.7e-3
        @test thickness(Filt) ≈ 0.15e-3
        @test thickness(Cover) ≈ 0.5e-3

        system = System([L1, L2, L3, Filt, Cover])

        # 0° beams
        beams = [
            Beam([0, -0.5e-3, -1.3e-3 / 2], [0, 1, 0], 0.5876e-6),
            Beam([0, -0.5e-3, 0], [0, 1, 0], 0.5876e-6),
            Beam([0, -0.5e-3, 1.3e-3 / 2], [0, 1, 0], 0.5876e-6)
        ]
        for beam in beams
            solve_system!(system, beam, r_max = 50)
            f_pos = last(beam.rays).pos + 0.12e-3 * last(beam.rays).dir

            # test if the beam is correctly focussed
            @test f_pos[3]≈0 atol=1e-7
        end

        # test rings
        @test 2*L1.shape.sdfs[4].hthickness ≈ 0.00060839 atol=1e-6
        @test 2*L2.shape.sdfs[4].hthickness ≈ 0.00057497 atol=1e-6
        @test 2*L3.shape.sdfs[4].hthickness ≈ 0.00048395 atol=1e-6
    end
end

@testset "Cylindrical Lenses" begin
    function working_distance(lens, offset_z)
        ray = Ray([0.0, -1.0, offset_z], [0.0, 1.0, 0])
        system = System([lens])
        beam = Beam(ray)
        solve_system!(system, beam)

        dist = -beam.rays[end].pos[3] / beam.rays[end].dir[3]
        α = asind(beam.rays[end].dir[3])
        wd = cosd(α) * dist

        return wd
    end
    @testset "Convex/cav cylinder lenses" begin
        # Thorlabs LJ1878L2, plano-convex
        r = 5.2e-3
        d = 10e-3
        h = 20e-3
        ct = 5.9e-3
        lens = Lens(
            CylindricalSurface(r, d, h),
            ct,
            n -> 1.517
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
        # test edge thickness
        @test thickness(BMO.shape(lens).sdfs[1])≈2.12e-3 atol=1e-4
        # test back focal length
        @test working_distance(lens, 0.05 * d / 2)≈6.1e-3 atol=1e-4

        # Thorlabs LK1900L1, plano-concave
        r = -13.1e-3
        d = 16e-3
        h = 18e-3
        ct = 2.0e-3
        lens = Lens(
            CylindricalSurface(r, d, h),
            ct,
            n -> 1.517
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
    end

    @testset "Convex/cav acylinder lenses" begin
        # Thorlabs AYL2520, plano-convex
        radius = 15.538e-3
        diameter = 25e-3
        height = 50e-3
        conic_constant = -1.0
        ct = 7.5e-3
        lens = Lens(
            BMO.AcylindricalSurface(
                radius,
                diameter,
                height,
                conic_constant,
                [0, 1.1926075e-5*(1e3)^3, -2.9323497e-9*(1e3)^5, -1.8718889e-11*(1e3)^7, -1.7009961e-14*(1e3)^9, 3.5481542e-17*(1e3)^11, 6.5241296e-20*(1e3)^13]
            ),
            ct,
            n -> 1.777
        )

        # test lens thickness
        @test thickness(lens) ≈ ct

        # test back focal length
        @test working_distance(lens, 0.05 * diameter / 2) ≈ 15.8e-3 atol=1e-4

        # Thorlabs AYL2520, plano-convex (inverted)
        radius = -15.538e-3
        diameter = 25e-3
        height = 50e-3
        conic_constant = -1.0
        ct = 7.5e-3
        lens = Lens(
            BMO.AcylindricalSurface(
                radius,
                diameter,
                height,
                conic_constant,
                [0, 1.1926075e-5*(1e3)^3, -2.9323497e-9*(1e3)^5, -1.8718889e-11*(1e3)^7, -1.7009961e-14*(1e3)^9, 3.5481542e-17*(1e3)^11, 6.5241296e-20*(1e3)^13]
            ),
            ct,
            n -> 1.777
        )

        # test lens thickness
        @test thickness(lens) ≈ ct
    end
end

@testset "Gaussian beamlet" begin
    @testset "Testing type definitions" begin
        @test isdefined(BMO, :GaussianBeamlet)
    end

    @testset "Testing analytical equations" begin
        λ = 500e-9
        w0 = 1e-3
        M2 = 1
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Test Rayleigh range and div. angle against Paschotta (https://www.rp-photonics.com/gaussian_beams.html)
        @test isapprox(zR, 6.28, atol = 1e-2)
        @test isapprox(BMO.beam_waist(zR, w0, zR), sqrt(2) * w0)
        @test isapprox(BMO.gouy_phase(zR, zR), -π / 4)
        @test isapprox(BMO.wavefront_curvature(zR, zR), 1 / (2 * zR))
        @test isapprox(BMO.divergence_angle(λ, w0, M2), 159e-6, atol = 1e-6)
    end

    @testset "Testing parameter correctness" begin
        # Gauss beam parameters
        y = -5:0.01:5       # m
        λ_1 = 500e-9        # m
        λ_2 = 1000e-9       # m
        P0 = 1              # W
        r = 0               # m
        w0_1 = 1e-3         # m
        w0_2 = 2e-3         # m
        M2_1 = 1e-3         # m
        M2_2 = 2e-3         # m
        E0_1 = BMO.electric_field(2 * P0 / (π * w0_1^2))
        E0_2 = BMO.electric_field(2 * P0 / (π * w0_2^2))
        gauss_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_1,
            w0_1,
            M2 = M2_1,
            P0 = P0)
        gauss_2 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ_2,
            w0_2,
            M2 = M2_2,
            P0 = P0)
        # Calculate analytical values
        zr_1 = BMO.rayleigh_range(λ_1, w0_1, M2_1)
        zr_2 = BMO.rayleigh_range(λ_2, w0_2, M2_2)
        wa_1 = BMO.beam_waist.(y, w0_1, zr_1)
        wa_2 = BMO.beam_waist.(y, w0_2, zr_2)
        Ra_1 = BMO.wavefront_curvature.(y, zr_1)
        Ra_2 = BMO.wavefront_curvature.(y, zr_2)
        ψa_1 = BMO.gouy_phase.(y, zr_1)
        ψa_2 = BMO.gouy_phase.(y, zr_2)
        Ea_1 = BMO.electric_field.(r, y, E0_1, w0_1, λ_1, M2_1)
        Ea_2 = BMO.electric_field.(r, y, E0_2, w0_2, λ_2, M2_2)
        # Calculate numerical values
        wn_1, Rn_1, ψn_1, w0n_1 = BMO.gauss_parameters(gauss_1, y)
        wn_2, Rn_2, ψn_2, w0n_2 = BMO.gauss_parameters(gauss_2, y)
        En_1 = [BMO.electric_field(gauss_1, r, yi) for yi in y]
        En_2 = [BMO.electric_field(gauss_2, r, yi) for yi in y]
        # Compare beam diameter within 0.1 nm
        @test all(isapprox.(wa_1, wn_1, atol = 1e-10))
        @test all(isapprox.(wa_2, wn_2, atol = 1e-10))
        # Compare wavefront curvature
        @test all(isapprox.(Ra_1, Rn_1, atol = 5e-9))
        @test all(isapprox.(Ra_2, Rn_2, atol = 5e-9))
        # Compare Gouy phase within 0.1 μrad
        @test all(isapprox.(ψa_1, ψn_1, atol = 1e-7))
        @test all(isapprox.(ψa_2, ψn_2, atol = 1e-7))
        # Compare calculated waist radius
        @test all(isapprox.(w0_1, w0n_1))
        @test all(isapprox.(w0_2, w0n_2))
        # Compare calculate electric field at r
        @test all(isapprox.(Ea_1, En_1, atol = 1e-8))
        @test all(isapprox.(Ea_2, En_2, atol = 1e-7))
        # Compare beam power with original value
        @test P0 ≈ BMO.optical_power(gauss_1)
        @test P0 ≈ BMO.optical_power(gauss_2)
    end

    @testset "Testing propagation correctness" begin
        # Analytical result using complex q factor
        q_ana(q0::Complex, M::Matrix) = (M[1] * q0 + M[3]) / (M[2] * q0 + M[4])
        R_ana(q::Complex) = real(1 / q)
        w_ana(q::Complex, λ, n = 1) = sqrt(-λ / (π * n * imag(1 / q)))
        propagate_ABCD(d) = [1 d; 0 1]
        lensmaker_ABCD(f) = [1 0; -1/f 1]
        # Beam parameters
        λ = 1000e-9
        w0 = 1e-3
        M2 = 1
        zr = BMO.rayleigh_range(λ, w0, M2)
        # Lens parameters
        R1 = 1
        R2 = 1
        lens_y_location = 0.1
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Stuff
        dy = 0.001
        ys = 0:dy:1.5
        w_analytical = Vector{Float64}(undef, length(ys))
        R_analytical = Vector{Float64}(undef, length(ys))
        # Propagate using ABCD formalism
        q0 = 0 + zr * im
        for i in 1:length(ys)
            w_analytical[i] = w_ana(q0, λ)
            R_analytical[i] = R_ana(q0)
            # catch first lens
            if i * dy == lens_y_location
                q0 = q_ana(q0, lensmaker_ABCD(f))
                continue
            end
            q0 = q_ana(q0, propagate_ABCD(dy))
        end

        # Numerical result
        tl = BMO.ThinLensSDF(R1, R2, 0.025)
        lens = Lens(tl, x -> nl)
        system = System(lens)
        translate3d!(lens, [0, lens_y_location, 0])
        # Create and solve beam, calculate beam parameters
        gauss = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
            λ,
            w0,
            support = [1, 0, 0],
            M2 = 1)
        solve_system!(system, gauss)
        w_numerical, R_numerical, ψ_numerical, w0_numerical = BMO.gauss_parameters(
            gauss,
            ys)
        # Compare beam radius to within 1 μm
        @test all(isapprox.(w_analytical, w_numerical, atol = 1e-6))
        # Compare if radius of curvature agreement to within 1 cm above 95% and any NaN
        temp = isapprox.(R_analytical, R_numerical, atol = 1e-2)
        @test sum(temp) / length(temp) > 0.95
        @test !any(isnan.(R_numerical))
        # Compare if Gouy phase zero at waists
        w0, i = findmin(w_analytical)
        @test isapprox(ψ_numerical[1], 0, atol = 1e-3)
        @test isapprox(ψ_numerical[i], 0, atol = 1e-3)
        # Compare calculated waist after lens
        @test isapprox(w0_numerical[i], w0, atol = 1e-7)

        @testset "Testing isparaxial and istilted" begin
            # Before lens rotation
            @test BMO.istilted(system, gauss) == false
            @test BMO.isparaxial(system, gauss) == true
            # Tilt lens, test again with 30° threshold for paraxial approx.
            zrotate3d!(lens, deg2rad(45))
            solve_system!(system, gauss)
            @test BMO.istilted(system, gauss) == true
            @test BMO.isparaxial(system, gauss, deg2rad(30)) == false
        end
    end
end

@testset "Detectors" begin
    @testset "Testing Spotdetector" begin
        # Set up tilted spot detection screens
        α = 45
        sd = Spotdetector(1.0)
        system = System([sd])
        translate3d!(sd, [0, 1, 0])
        zrotate3d!(sd, deg2rad(α))
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        # Trace beams in x-y-plane
        xs = LinRange(-0.25, 0.25, 10)
        for x in xs
            BMO.position!(first(beam.rays), Point3{Float64}(x, 0, 0))
            solve_system!(system, beam)
            # compare ray intersection to stored data
            data = last(sd.data)
            ray = last(beam.rays)
            pos_ray = position(ray) +
                      length(BMO.intersection(ray)) * BMO.direction(ray)
            pos_dta = position(sd) +
                      BMO.orientation(sd)[:, 1] * data[1]
            @test pos_ray ≈ pos_dta
        end
        # Test reset function
        empty!(sd)
        @test isempty(sd.data)
    end
end

@testset "Interference" begin
    @testset "Pre-Beamsplitter tests with seperate beams" begin
        # Gauss beam parameters (selected for ring fringes)
        w0 = 0.01e-3
        λ = 1000e-9
        M2 = 1
        P0 = 1e-3
        I0 = 2 * P0 / (π * w0^2)
        E0 = BMO.electric_field(I0)
        zR = BMO.rayleigh_range(λ, w0, M2)
        # Detector parameters
        z = 0.1     # distance to detector
        l = 1e-2    # detector size
        n = 1000    # detector grid resolution
        # Lens parameters
        R1 = R2 = d = 0.01
        nl = 1.5
        f = BMO.lensmakers_eq(R1, -R2, nl)
        # Raytracing system (for all tests)
        pd_l = Photodetector(l, n)
        pd_s = Photodetector(l / 10, n ÷ 10)
        ln = ThinLens(R1, R2, d, nl)
        translate3d!(pd_l, [0, z, 0])
        translate3d!(pd_s, [0, z, 0])
        translate3d!(ln, [0, z - f - thickness(ln.shape) / 2, 0])

        @testset "Testing fringe pattern" begin
            system = System(pd_l)
            Δz = 5e-3   # arm length difference
            # Analytic solution
            xs = ys = LinRange(-l / 2, l / 2, n)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(r, z, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(r, z + Δz, E0, w0, λ, M2)
                end
            end
            # Numerical solution
            empty!(pd_l)
            g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            g_2 = GaussianBeamlet([0.0, -Δz, 0], [0.0, 1, 0],
                λ,
                w0,
                M2 = M2,
                P0 = P0)
            solve_system!(system, g_1)
            solve_system!(system, g_2)

            # Compare solutions
            I_analytical = BMO.intensity.(screen)
            I_numerical = BMO.intensity.(pd_l.field)
            Pt = BMO.optical_power(pd_l)
            @test all(isapprox.(I_analytical, I_numerical, atol = 2e-1))
            @test isapprox(Pt, 2 * P0, atol = 3e-5)
        end

        @testset "Testing λ phase shift" begin
            system = System([pd_s, ln])
            # Numerical solution
            Δz = LinRange(0, λ, 50)
            Pt_numerical = zeros(length(Δz))
            for (i, z_i) in enumerate(Δz)
                empty!(pd_s)
                g_1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                g_2 = GaussianBeamlet([0.0, z_i, 0], [0.0, 1, 0],
                    λ,
                    w0,
                    M2 = M2,
                    P0 = P0)
                solve_system!(system, g_1)
                solve_system!(system, g_2)
                Pt_numerical[i] = BMO.optical_power(pd_s)
                # Test length/opl function
                @test length(g_1) == z
                @test length(g_1) < BMO.optical_path_length(g_1)
                @test length(g_2) == length(g_1) - z_i
            end
            # Analytical solution (cosine over Δz), ref. power is 4*P0 since beamsplitter is missing
            Pt_analytical = 4 * P0 * [(cos(2π * z / (maximum(Δz))) + 1) / 2 for z in Δz]

            # Compare detectors (this also tests correct behavior when focussing the beam)
            @test all(isapprox.(Pt_numerical, Pt_analytical, atol = 1e-4))
        end
    end

    @testset "Michelson Interferometer" begin
        # setup Michelson Interferometer
        l_0 = 0.1
        pd_size = BMO.inch / 5
        pd_resolution = 100
        m1 = SquarePlanoMirror2D(BMO.inch)
        m2 = SquarePlanoMirror2D(BMO.inch)
        bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
        pd = Photodetector(pd_size, pd_resolution)
        translate3d!(m1, [l_0, 0, 0])
        translate3d!(m2, [0, l_0, 0])
        translate3d!(pd, [-l_0, 0, 0])
        zrotate3d!(bs, deg2rad(45))
        zrotate3d!(m1, deg2rad(90))
        zrotate3d!(pd, deg2rad(90))

        system = System([m1, m2, bs, pd])

        # Test correct values for reflectivity/transmission
        @test isvalid(bs)

        @testset "Equal armlength MI - integrated power" begin
            # setup 635 nm laser with 0.1 mm waist for fast divergence
            λ = 635e-9
            P_0 = 5e-3
            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, 1e-4, P0 = P_0)

            # Shift mirror #2 by -λ to +λ
            lambdas = LinRange(-λ, λ, 200)

            path_length_numerical = zeros(length(lambdas))
            optical_pwr_numerical = zeros(length(lambdas))

            for (i, lambda) in enumerate(lambdas)
                translate_to3d!(m2, [0, l_0, 0] + [0, lambda, 0])
                empty!(pd)
                solve_system!(system, beam)

                # Moving mirror path length
                path_length_numerical[i] = length(beam.children[1].children[2])
                optical_pwr_numerical[i] = BMO.optical_power(pd)
            end

            path_length_analytical = @. 2 * lambdas + 4l_0
            optical_pwr_analytical = @. P_0 * (1 / 2 * cos(2π * (2lambdas / λ) + π) + 1 / 2)

            # Compare correct PD signal and λ shift in moving arm
            @test all(isapprox.(optical_pwr_analytical, optical_pwr_numerical, atol = 5e-6))
            @test all(isapprox.(path_length_analytical, path_length_numerical))
        end

        @testset "Unequal armlength MI - electrical field" begin
            λ = 635e-9
            w0 = 1e-4
            P0 = 1e-3
            M2 = 1
            I0 = 2 * P0 / (π * w0^2)
            E0 = BMO.electric_field(I0) * 1 / sqrt(2)^2
            zR = BMO.rayleigh_range(λ, w0, M2)

            beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], λ, w0, P0 = P0, M2 = M2)

            # arm length diff
            Δl = 1 * l_0
            translate_to3d!(m2, [0, l_0 + Δl, 0])

            # numerical solution
            empty!(pd)
            solve_system!(system, beam)

            # analytical solution
            short_arm = 4l_0
            long_arm = short_arm + 2Δl
            xs = ys = LinRange(-pd_size / 2, pd_size / 2, pd_resolution)
            screen = zeros(ComplexF64, length(xs), length(ys))
            for (j, y) in enumerate(ys)
                for (i, x) in enumerate(xs)
                    r = sqrt(x^2 + y^2)
                    screen[i, j] += BMO.electric_field(
                        r, short_arm, E0, w0, λ, M2)
                    screen[i, j] += BMO.electric_field(
                        r, long_arm, E0, w0, λ, M2) * exp(im * pi)
                end
            end

            Re_analytical = real.(screen)
            Im_analytical = imag.(screen)

            Re_numerical = real(pd.field)
            Im_numerical = imag(pd.field)

            # Compare solutions, units V/m
            @test all(isapprox.(Re_analytical, Re_numerical, atol = 5e-2))
            @test all(isapprox.(Im_analytical, Im_numerical, atol = 5e-2))
        end
    end

    @testset "Testing power conservation" begin
        # variables
        P0 = 0.5 # W
        l0 = 0.1 # m
        w0 = 0.5e-3
        λ = 1064e-9

        bs = ThinBeamsplitter(10e-3)
        pd_1 = Photodetector(10e-3, 100)
        pd_2 = Photodetector(10e-3, 100)

        zrotate3d!(bs, deg2rad(45))
        translate3d!(pd_1, [0, l0, 0])
        zrotate3d!(pd_1, deg2rad(180))

        translate3d!(pd_2, [l0, 0, 0])
        zrotate3d!(pd_2, deg2rad(90))

        # add BS and PD orientation error
        zrotate3d!(bs, deg2rad(0.017))
        zrotate3d!(pd_1, deg2rad(10))
        xrotate3d!(pd_1, deg2rad(15))

        # define system and beams -> solve
        system = System([bs, pd_1, pd_2])

        phis = LinRange(0, 2pi, 25)
        p1 = similar(phis)
        p2 = similar(phis)

        l1 = GaussianBeamlet([0, -l0, 0], [0, 1.0, 0], λ, w0; P0)
        l2 = GaussianBeamlet([-l0, 0, 0], [1.0, 0, 0], λ, w0; P0)

        E0_buffer = l1.E0

        for (i, phi) in enumerate(phis)
            # Iterate over relative phase shifts, use retracing
            l1.E0 = E0_buffer * exp(im * phi)
            empty!(pd_1)
            empty!(pd_2)
            solve_system!(system, l1)
            solve_system!(system, l2)
            p1[i] = BMO.optical_power(pd_1)
            p2[i] = BMO.optical_power(pd_2)
            # Test power conservation
            @test p1[i] + p2[i] - 2P0 < 1e-4 # W
        end
    end
end

@testset "Polarized rays" begin
    @testset "Polarization transforms" begin
        # Reflection matrix and lin. x-pol
        J = [-1 0 0; 0 1 0; 0 0 1]
        E0 = [1, 0, 0]

        @testset "90° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [1, 0, 0]
            @test BMO._calculate_global_E0(in_dir, out_dir, J, E0) ≈ [0, 0, -1]
        end

        @testset "0° reflection" begin
            in_dir = [0, 0, 1]
            out_dir = [0, 0, -1]
            @test BMO._calculate_global_E0(in_dir, out_dir, J, E0) ≈ [-1, 0, 0]
        end
    end

    @testset "Mirror reflections" begin
        # Setup system as in https://opg.optica.org/ao/fulltext.cfm?uri=ao-50-18-2855&id=218813
        m1 = SquarePlanoMirror2D(1.0)
        m2 = SquarePlanoMirror2D(1.0)
        m3 = SquarePlanoMirror2D(1.0)
        translate3d!(m2, [2, 0, 0])
        translate3d!(m3, [2, 2, 0])
        zrotate3d!(m1, deg2rad(-90))
        yrotate3d!(m1, deg2rad(45))
        zrotate3d!(m2, deg2rad(45))
        xrotate3d!(m3, deg2rad(135))

        system = StaticSystem([m1, m2, m3])

        I0_1 = 1
        I0_2 = 5
        lin_x_pol = [I0_1, 0, 0]
        lin_y_pol = [0, I0_2, 0]

        # Beam of polarized rays
        ray = PolarizedRay([0.0, 0, -2], [0, 0, 1], 1000e-9, lin_x_pol)
        beam = Beam(ray)

        @testset "x-Polarization" begin
            BMO.polarization!(ray, lin_x_pol)
            # test tracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_x_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, 0, -I0_1]
            @test BMO.polarization(beam.rays[3]) ≈ [0, 0, I0_1]
            @test BMO.polarization(beam.rays[4]) ≈ [0, -I0_1, 0]
            @test length(beam) == 6.0
        end

        @testset "y-Polarization" begin
            BMO.polarization!(ray, lin_y_pol)
            translate3d!(m3, [0, 2, 0])
            # test retracing
            solve_system!(system, beam)
            @test BMO.polarization(beam.rays[1]) ≈ lin_y_pol
            @test BMO.polarization(beam.rays[2]) ≈ [0, -I0_2, 0]
            @test BMO.polarization(beam.rays[3]) ≈ [I0_2, 0, 0]
            @test BMO.polarization(beam.rays[4]) ≈ [-I0_2, 0, 0]
            @test length(beam) == 8.0
        end
    end

    @testset "Brewster windows" begin
        brewster_angle(n) = atan(n)
        # Testcase based on 5 successive Brewster windows
        n = 1.5
        θb = brewster_angle(n)
        d = 0.1
        # Calculate transmission efficiency
        rs, rp, ts, tp = BMO.fresnel_coefficients(θb, n)
        Ts = 1 - abs2(rs)
        Tp = 1 - abs2(rp)
        # Setup testcase
        s1 = BMO.CuboidMesh(1.0, d, 1.0)
        s2 = BMO.CuboidMesh(1.0, d, 1.0)
        s3 = BMO.CuboidMesh(1.0, d, 1.0)
        s4 = BMO.CuboidMesh(1.0, d, 1.0)
        s5 = BMO.CuboidMesh(1.0, d, 1.0)
        l1 = Lens(s1, x -> n)
        l2 = Lens(s2, x -> n)
        l3 = Lens(s3, x -> n)
        l4 = Lens(s4, x -> n)
        l5 = Lens(s5, x -> n)
        translate3d!.([l1, l2, l3, l4, l5], Ref([-0.5, -d / 2, -0.5]))
        BMO.set_new_origin3d!.(BMO.shape.([l1, l2, l3, l4, l5]))
        translate3d!(l2, [0, 0.5, -1d / 2])
        translate3d!(l3, [0, 1.0, -2d / 2])
        translate3d!(l4, [0, 1.5, -3d / 2])
        translate3d!(l5, [0, 2.0, -4d / 2])
        xrotate3d!.([l1, l2, l3, l4, l5], -θb)
        # Solve system of s- and p-polarized beams
        system = StaticSystem([l1, l2, l3, l4, l5])
        x_pol_ray = PolarizedRay(
            [-0.1, -1, 0], [0, 1.0, 0], 1000e-9, [BMO.electric_field(1), 0, 0])
        z_pol_ray = PolarizedRay(
            [+0.1, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        s_beam = Beam(x_pol_ray)
        p_beam = Beam(z_pol_ray)
        solve_system!(system, s_beam)
        solve_system!(system, p_beam)
        # Since system is non-focussing, calculate pseudo-intensity
        pseudo_Is = abs2(BMO.polarization(last(BMO.rays(s_beam)))[1]) /
                    (2 * BMO.Z_vacuum)
        pseudo_Ip = abs2(BMO.polarization(last(BMO.rays(p_beam)))[3]) /
                    (2 * BMO.Z_vacuum)
        # Test against m interfaces
        m = length(system.objects) * 2
        @test pseudo_Is ≈ Ts^m
        @test pseudo_Ip ≈ Tp^m
    end

    @testset "Fresnel rhomb" begin
        # Create Fresnel rhomb with n=1.5 and θ=53.3° for quarter-wave plate effect
        n = 1.5
        s1 = BMO.CuboidMesh(0.5, 1.25, 0.5, deg2rad(53.3))
        l1 = Lens(s1, x -> n)
        translate3d!(l1, [-0.25, 0, -0.25])
        BMO.set_new_origin3d!(s1)
        # Rotate prism to obtain 45° beam input polarization
        yrotate3d!(l1, deg2rad(135))
        # Solve system
        system = StaticSystem([l1])
        ray = PolarizedRay(
            [0, -1, 0], [0, 1.0, 0], 1000e-9, [0, 0, BMO.electric_field(1)])
        beam = Beam(ray)
        solve_system!(system, beam)
        # Assumes propagation along the y-axis after rhomb, calculate polarization state
        Ex = getindex.(BMO.polarization.(beam.rays), 1)
        Ey = getindex.(BMO.polarization.(beam.rays), 2)
        Ez = getindex.(BMO.polarization.(beam.rays), 3)
        # Test for circular polarization and Ey error
        phi = angle(last(Ez)) - angle(last(Ex))
        @test phi ≈ π / 2
        @test abs(last(Ey)) < 2e-14
    end

    @testset "Mach-Zehnder Interferometer" begin
        # setup MZI
        m1 = SquarePlanoMirror2D(BMO.inch)
        m2 = SquarePlanoMirror2D(BMO.inch)
        b1 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
        b2 = ThinBeamsplitter(BMO.inch, reflectance = 0.5)

        system = StaticSystem([m1, m2, b1, b2])

        translate3d!(b1, [0 * BMO.inch, 0 * BMO.inch, 0])
        translate3d!(b2, [2 * BMO.inch, 2 * BMO.inch, 0])
        translate3d!(m1, [0 * BMO.inch, 2 * BMO.inch, 0])
        translate3d!(m2, [2 * BMO.inch, 0 * BMO.inch, 0])

        # Rotate with consideration to mirror/bs normal
        zrotate3d!(b1, deg2rad(360 - 135))
        zrotate3d!(b2, deg2rad(45))
        zrotate3d!(m1, deg2rad(360 - 135))
        zrotate3d!(m2, deg2rad(45))

        ray = PolarizedRay([0, -0.1, 0], [0.0, 1.0, 0], 1000e-9, [0, 0, 1])
        beam = Beam(ray)

        @testset "z-polarized ray along y-axis" begin
            # Solve with z-polarized ray along y-axis
            BMO.polarization!(ray, [0, 0, 1])
            solve_system!(system, beam)

            # Extract E0s: t - transmitted, r - reflected
            t = beam.children[1].rays[1].E0
            r = beam.children[2].rays[1].E0
            tr = beam.children[1].rays[2].E0
            rr = beam.children[2].rays[2].E0
            trt = beam.children[1].children[1].rays[1].E0
            trr = beam.children[1].children[2].rays[1].E0
            rrt = beam.children[2].children[1].rays[1].E0
            rrr = beam.children[2].children[2].rays[1].E0

            # Test phase flips
            @test t[3] ≈ sqrt(2) / 2
            @test r[3] ≈ -sqrt(2) / 2
            @test tr[3] ≈ -sqrt(2) / 2
            @test rr[3] ≈ sqrt(2) / 2
            @test trt ≈ rrr
            @test trr ≈ rrt
        end

        @testset "x-polarized ray along y-axis" begin
            # Test num. of leaves before retracing
            @test length(collect(Leaves(beam))) == 4
            # Retrace with z-polarized ray along y-axis
            BMO.polarization!(ray, [1, 0, 0])
            solve_system!(system, beam)

            # Extract E0s: t - transmitted, r - reflected
            t = beam.children[1].rays[1].E0
            r = beam.children[2].rays[1].E0
            tr = beam.children[1].rays[2].E0
            rr = beam.children[2].rays[2].E0
            trt = beam.children[1].children[1].rays[1].E0
            trr = beam.children[1].children[2].rays[1].E0
            rrt = beam.children[2].children[1].rays[1].E0
            rrr = beam.children[2].children[2].rays[1].E0

            # Test phase flips
            @test t[1] ≈ sqrt(2) / 2
            @test r[2] ≈ -sqrt(2) / 2
            @test tr ≈ r
            @test rr ≈ t
            @test trt ≈ rrr
            @test trr ≈ rrt
        end
    end
end

@testset "Beamsplitters" begin
    mm = 1e-3
    N0 = 1.5
    @testset "Testing RectangularPlateBeamsplitter with Beam" begin
        # Init splitter
        N0 = 1.5
        mm = 1e-3
        pbs = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, n -> N0)
        system = System([pbs])
        beam = Beam([0, -50mm, 0], [0, 1, 0], 1e-6)
        # Trace normally
        zrotate3d!(pbs, deg2rad(45))
        solve_system!(system, beam)

        @testset "Test pos/dir" begin
            @test position(pbs) == zeros(3)
            @test BMO.orientation(pbs) ≈ BMO.orientation(pbs.substrate)
        end

        @testset "Test children after tracing" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 1
            @test length(t) == 2
            @test length(r) == 1
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== 1)
            @test all(BMO.refractive_index.(t) .== [N0, 1])
            @test all(BMO.refractive_index.(r) .== 1)
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(first(r)) ≈ [1, 0, 0]
        end

        # Retrace backside
        zrotate3d!(pbs, π)
        solve_system!(system, beam)

        @testset "Test children after retracing" begin
            p = beam.rays
            t = beam.children[1].rays
            r = beam.children[2].rays
            # no of rays
            @test length(p) == 2
            @test length(t) == 1
            @test length(r) == 2
            # correct ref. index
            @test all(BMO.refractive_index.(p) .== [1, N0])
            @test all(BMO.refractive_index.(t) .== 1)
            @test all(BMO.refractive_index.(r) .== [N0, 1])
            # correct dir
            @test BMO.direction(first(p)) ≈ BMO.direction(last(t))
            @test BMO.direction(last(r)) ≈ [1, 0, 0]
        end
    end

    @testset "Testing CubeBeamsplitter with Beam" begin
        # Init splitter
        cbs = CubeBeamsplitter(25e-3, n -> N0)
        translate3d!(cbs, [0, 50mm, 0])
        system = System([cbs])
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)

        @testset "Initial CBS tracing" begin
            # Trace normally
            solve_system!(system, beam)
            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end

        @testset "Retrace after 45° CBS rotation" begin
            # Retrace
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test BMO.direction(last(t)) == BMO.direction(first(p))
            @test BMO.direction(last(t)) == [0, 1, 0]
        end

        @testset "Retrace CBS backside" begin
            # Retrace backside
            zrotate3d!(cbs, π / 2)
            solve_system!(system, beam)

            # Test correct ray length, ref. indices, dirs
            p = BMO.rays(beam)
            t = BMO.rays(beam.children[1])
            r = BMO.rays(beam.children[2])

            @test length(p) == 2
            @test length(r) == 2
            @test length(t) == 2
            @test BMO.refractive_index.(p) == [1, N0]
            @test BMO.refractive_index.(t) == [N0, 1]
            @test BMO.refractive_index.(r) == [N0, 1]
            @test BMO.direction(last(t)) ≈ BMO.direction(first(p))
            @test BMO.direction(last(r)) ≈ [-1, 0, 0]
        end
    end
end

@testset "Dummy objects" begin
    # Setup dummy cube and test beam
    cube_shape = BMO.CubeMesh(1)
    translate3d!(cube_shape, -[0.5, 0, 0.5])
    translate3d!(cube_shape, [0, 5, 0])
    @testset "IntersectableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        intersectable = IntersectableObject(cube_shape)
        system = System([intersectable])
        solve_system!(system, beam)
        # Test nothing interaction
        @test length(BMO.rays(beam)) == 1
        @test BMO.object(BMO.intersection(last(BMO.rays(beam)))) == intersectable
        @test BMO.shape(BMO.intersection(last(BMO.rays(beam)))) == cube_shape
        @test isnothing(BMO.interact3d(system, intersectable, beam, first(BMO.rays(beam))))
    end

    @testset "NonInteractableObject" begin
        beam = Beam([0, 0, 0], [0, 1, 0], 1e-6)
        noninteract = NonInteractableObject(cube_shape)
        system = System([noninteract])
        solve_system!(system, beam)
        # Test nothing interaction and intersection
        @test length(BMO.rays(beam)) == 1
        @test isnothing(BMO.intersection(last(BMO.rays(beam))))
        @test isnothing(BMO.interact3d(system, noninteract, beam, first(BMO.rays(beam))))
    end
end

@testset "Render" begin
    axis = nothing
    cube = BMO.CubeMesh(1)
    @test_throws BMO.MissingBackendError render!(axis, cube)
end

@testset "Double Gauss lens" begin
    """Tests hit distance to detector origin"""
    function test_coma(pd::Spotdetector; atol=7e-5)
        for hit in pd.data
            if norm(hit) > atol
                return false
            end
        end
        return true
    end
    
    # Based on https://www.pencilofrays.com/double-gauss-sonnar-comparison/
    l1 = SphericalLens(48.88e-3, 182.96e-3, 8.89e-3, 52.3e-3, λ -> 1.62286)
    l23 = SphericalDoubletLens(36.92e-3, Inf, 23.06e-3, 15.11e-3, 2.31e-3,
    45.11e-3, λ -> 1.58565, λ -> 1.67764)
    l45 = SphericalDoubletLens(-23.91e-3, Inf, -36.92e-3, 1.92e-3, 7.77e-3,
    40.01e-3, λ -> 1.57046, λ -> 1.64128)
    l6 = SphericalLens(1063.24e-3, -48.88e-3, 6.73e-3, 45.11e-3, λ -> 1.62286)
    # Calculate translation distances
    l_23 = thickness(l1) + 0.38e-3
    l_45 = l_23 + thickness(l23) + 9.14e-3 + 13.36e-3
    l_6 = l_45 + thickness(l45) + 0.38e-3
    # Corresponds to back focal length of f=59.21 mm on y-axis from link above + "error" δf
    δf = 7e-4
    f_z = l_6 + thickness(l6.shape) + 58.21e-3 + δf
    translate3d!(l23, [0, l_23, 0])
    translate3d!(l45, [0, l_45, 0])
    translate3d!(l6, [0, l_6, 0])
    # Create and move group - this tests a bunch of kinematic correctness
    double_gauss = ObjectGroup([l1, l23, l45, l6])
    
    # Spawn Spotdetector, move to "focus"
    detector = Spotdetector(5e-3)
    translate3d!(detector, [0,f_z,0])    
    # Setup test, rotate and translate to test correct kinematics
    test_setup = ObjectGroup([double_gauss, detector])
    translate3d!(test_setup, [0.05, 0.05, 0.05])
    xrotate3d!(test_setup, deg2rad(60))
    zrotate3d!(test_setup, deg2rad(45))    
    system = System([test_setup])
    
    @testset "Test with collimated source" begin    
        ## Test against back focal length as per source above
        dir = BMO.orientation(double_gauss)[:, 2] # rotated collimated ray direction
        pos = position(l1) - 0.05 * dir # rotated collimated ray position    
        source = CollimatedSource(pos, dir, 0.04, 486.0e-9, num_rays=1000, num_rings=10)
        solve_system!(system, source)    
        @test test_coma(detector, atol=2e-5)
    end
    
    # Reset test scenario back to origin, y-axis alignment
    reset_rotation3d!(test_setup)
    reset_translation3d!(test_setup)
    translate_to3d!(detector, [0, 0.147, 0])
    
    @testset "Test with point source (wide)" begin
        source = PointSource([0, -0.5, 0], [0, 1 ,0], deg2rad(2), 486.0e-9, num_rays=1000, num_rings=10)
        empty!(detector)
        solve_system!(system, source)
        @test test_coma(detector, atol=6e-5)
    end
    
    @testset "Test with point source (narrow)" begin
        # Tests regression for https://github.com/StackEnjoyer/BeamletOptics.jl/issues/11
        source = PointSource([0, -0.5, 0], [0, 1 ,0], 5e-5, 486.0e-9, num_rays=1000, num_rings=10)
        empty!(detector)
        solve_system!(system, source)
        @test test_coma(detector, atol=2e-7)
    end
end

@testset "Bug fixes" begin
    @testset "Issue#14" begin
        pd_res = 1000
        pd = BMO.Photodetector(10e-3, pd_res)
        BMO.zrotate3d!(pd, deg2rad(90))
        BMO.translate3d!(pd, [0.46, 0, 0])
        # Setup beam
        y_0 = 0.2
        beam = BMO.GaussianBeamlet(
            [0, y_0, 0], [0.46, -y_0, 0], 532e-9, 2.5e-3, P0 = 10e-3)
        # Solve system
        system = BMO.System([pd])
        empty!(pd)
        BMO.solve_system!(system, beam)

        @test BMO.optical_power(pd)≈10e-3 atol=1e-5
    end
end

@testset "Aqua" begin
    using Aqua
    Aqua.test_all(BMO)
end