module TestUtils

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Utilities" begin
    @testset "Testing normal3d" begin
        @testset "normal3d(v) with Array" begin
            v = [1.0, 1, 1]
            k = BMO.normal3d(v)
            @test dot(v, k)≈0 atol=2e-14
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
 
    @testset "Testing isparallel3d and isorthogonal3d" begin
        v1 = [1,0,0]
        v2 = Point3(1,0,eps())
        v3 = BMO.normal3d(v1)
        # test parallel
        @test BMO.isparallel3d(v1, v2)
        @test !BMO.isparallel3d(v1, v3)
        # test orthogonal
        @test !BMO.isorthogonal3d(v1, v2)
        @test BMO.isorthogonal3d(v1, v3)
    end

    @testset "Testing rotate3d for clockwise dir. and conservation of length" begin
        Rot = BMO.rotate3d([0, 0, 1], π / 2)
        @test isapprox(Rot * [1, 0, 0], [0, 1, 0])
        @test_throws ArgumentError BMO.rotate3d([0, 0, 1], Inf)
        @test_throws ArgumentError BMO.rotate3d([0, 0, 1], NaN)
    end

    @testset "Testing align3d for rotation and conservation of length" begin
        # Start vector must have unit length!
        start = [1, 0, 0]
        # Test parallel case
        target = [1, 0, 0]
        T = BMO.align3d(start, target)
        @test T * start ≈ target
        # Test parallel opposite cases (X, Y, Z axes and arbitrary 3D direction)
        for s_vec in ([1, 0, 0], [0, 1, 0], [0, 0, 1], normalize([1, 2, 3]))
            T_opp = BMO.align3d(s_vec, -s_vec)
            @test T_opp * s_vec ≈ -s_vec
            @test det(T_opp) ≈ 1
        end
        # Test norm and 45° rotation
        target = [1.0, 1.0, 0.0]
        T = BMO.align3d(start, target)
        @test T * start ≈ normalize(target)
    end

    @testset "Testing angle3d for resulting angle" begin
        a = BMO.angle3d([1, 0, 0], [0, 0, 1])
        @test isapprox(a, π / 2)
        # test signed angles
        b1 = BMO.angle3d([1,0,0], [0,1,0], [0,0,1])
        b2 = BMO.angle3d([1,0,0], [0,-1,0], [0,0,1])
        b3 = BMO.angle3d([1,1,0], [1,0,1], -[0,0,1])
        @test b1 ≈ -π/2
        @test b2 ≈ π/2
        @test b3 ≈ -deg2rad(60)
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
        ref_sellm = SellmeierEquation(0.6961663, 0.4079426, 0.8974794,
                                      0.0684043^2, 0.1162414^2, 9.896161^2)

        @testset "DiscreteRefractiveIndex" begin
            @test isdefined(BMO, :DiscreteRefractiveIndex)
            @test isa(ref_index, DiscreteRefractiveIndex{T2})
            @test ref_index(lambdas[2]) == indices[2]
            @test_throws KeyError ref_index(lambdas[1] + 1e-9)
            # Test constructor
            @test_throws ArgumentError DiscreteRefractiveIndex([1], [1, 2])
        end

        @testset "SellmeierEquation" begin
            # Known reference index for fused silica at 500 nm
            @test isapprox(ref_sellm(500e-9), 1.4623, atol=3e-5)
            # Check that wavelength in μm gives same result as m (via internal conversion)
            λ_m = 1.0e-6
            λ_um = 1.0
            @test isapprox(ref_sellm(λ_m), ref_sellm(λ_um * 1e-6), atol=1e-12)
            # Check monotonic decrease of n with λ in the visible range
            n_vis = [ref_sellm(λ) for λ in (400:50:700) .* 1e-9]
            @test all(diff(n_vis) .< 0)
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
            @test isnothing(BMO.test_refractive_index_function(ref_sellm))
        end
    end

    @testset "Testing numerical aperture util." begin
        n0 = 1.5
        theta = deg2rad(45)
        BMO.numerical_aperture(theta, n0) == n0 * sin(theta)
    end

    @testset "Testing bisection finder" begin
        f(x) = x^2 + 1
        g(x) = x^3 + x^2
        @test_throws ErrorException BMO.find_zero_bisection(f, -1, 1)
        @test_throws ErrorException BMO.find_zero_bisection(g, -2, 1; tol=5e-16, max_iter=10)
    end

    @testset "Testing interferometric visibility util." begin
        Iphi(phi, V) = V*sin(phi) + 1 
        V = 0.85
        pwr = [Iphi(phi, V) for phi = LinRange(0, 2pi, 1000)]
        @test BMO.visibility(pwr) ≈ V atol=2e-6
    end

    @testset "Testing polarization state utils." begin
        lin1 = [1, 0, 0]
        lin2 = [-1, 1, 0]
        circ = normalize([1,1*exp(im*pi/2),0])
        ellp = [1, exp(im*pi/2) * 2, 0]
        
        @test BMO.islinear(lin1)
        @test BMO.islinear(lin2)
        @test !BMO.islinear(circ)
        @test !BMO.islinear(ellp)
        
        @test !BMO.iscircular(lin1)
        @test !BMO.iscircular(lin2)
        @test BMO.iscircular(circ)
        @test !BMO.iscircular(ellp)
        
        @test !BMO.iselliptical(lin1)
        @test !BMO.iselliptical(lin2)
        @test !BMO.iselliptical(circ)
        @test BMO.iselliptical(ellp)
    end

    @testset "Testing ellipse parameterization" begin
        ts = LinRange(0, 2pi, 1000)
        p0 = Point2(1,1)
        r1 = [0,1]
        # 45° deg inclination
        r2 = [1,1.]
        pts = BMO.ellipse(ts, p0, r1, r2)
        x = getindex.(pts, 1)
        z = getindex.(pts, 2)
        # test min/max values of 2D bounding box
        @test isapprox(maximum(x), 2, atol=1e-5)
        @test isapprox(minimum(x), 0, atol=1e-5)
        @test isapprox(maximum(z), 1+sqrt(2), atol=1e-5)
        @test isapprox(minimum(z), 1-sqrt(2), atol=1e-5)
    end
end

end # MODULE