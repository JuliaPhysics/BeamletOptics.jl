module TestBeamGroups

using BeamletOptics
using Test
using LinearAlgebra
using GeometryBasics

const BMO = BeamletOptics

@testset "Beam groups" begin
    @testset "AbstractBeamGroup definitions" begin
        @test isdefined(BMO, :AbstractBeamGroup)

        mutable struct BeamTestGroup{T} <: BMO.AbstractBeamGroup{T, Ray{T}}
            central_beam::BMO.Beam{T, Ray{T}}
            center::Point3{T}
            orientation::BMO.SMatrix{3, 3, T, 9}
        end

        struct BeamTestSystem <: BMO.AbstractSystem end
        BMO.objects(::BeamTestSystem) = [nothing]
        BMO.intersect3d(::Nothing, ::BMO.AbstractRay) = nothing
        BMO.interact3d(::Nothing) = nothing

        BMO.beams(tg::BeamTestGroup) = [tg.central_beam]

        pos = [0,0,0]
        dir = [0,1,0]
        lambda = 1064e-9

        tg = BeamTestGroup{Float64}(BMO.Beam(pos, dir, lambda), pos, BMO.SMatrix{3, 3, Float64, 9}(I))

        @test position(tg) == pos
        @test BMO.direction(tg) == dir
        @test orientation(tg) == I
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

        @testset "Testing basis kwarg" begin
            dir_n = normalize(dir)
            # a basis is only defined up to its component normal to dir
            b0 = BMO.normal3d(dir_n)
            rotated = PointSource(pos, dir, alpha, lambda; num_rays, num_rings,
                basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)

            dirs_default = BMO.direction.(first.(BMO.rays.(BMO.beams(source))))
            dirs_rotated = BMO.direction.(first.(BMO.rays.(BMO.beams(rotated))))

            # rotating the basis spins the fan about its own axis: the polar angle
            # distribution is invariant, the individual ray directions are not
            angles_default = sort(round.(BMO.angle3d.(Ref(dir), dirs_default), digits=11))
            angles_rotated = sort(round.(BMO.angle3d.(Ref(dir), dirs_rotated), digits=11))
            @test angles_default ≈ angles_rotated
            @test !all(dirs_default .≈ dirs_rotated)

            # the default is reproducible, and passing the default basis reproduces it
            @test all(dirs_default .≈ BMO.direction.(first.(BMO.rays.(BMO.beams(
                PointSource(pos, dir, alpha, lambda; num_rays, num_rings))))))
            @test all(dirs_default .≈ BMO.direction.(first.(BMO.rays.(BMO.beams(
                PointSource(pos, dir, alpha, lambda; num_rays, num_rings, basis = b0))))))

            # a basis component along dir is projected out, so it changes nothing
            @test all(dirs_default .≈ BMO.direction.(first.(BMO.rays.(BMO.beams(
                PointSource(pos, dir, alpha, lambda; num_rays, num_rings,
                    basis = b0 + 5 * dir_n))))))
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException PointSource(pos, dir, 1.1*π, lambda; num_rays, num_rings)
            @test_throws ErrorException PointSource(pos, dir, alpha, lambda; num_rays=100, num_rings=10)
            # basis parallel to dir has no component in the sampling plane
            @test_throws ErrorException PointSource(pos, dir, alpha, lambda; num_rays, num_rings,
                basis = dir)
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

        @testset "Testing basis kwarg" begin
            b0 = BMO.normal3d(dir)
            rotated = CollimatedSource(pos, dir, diameter; num_rays, num_rings,
                basis = BMO.rotate3d(dir, deg2rad(30)) * b0)

            pos_default = position.(first.(BMO.rays.(BMO.beams(source))))
            pos_rotated = position.(first.(BMO.rays.(BMO.beams(rotated))))

            # rotating the basis spins the pupil pattern about its own axis: the radial
            # distribution is invariant, the individual ray positions are not
            radii_default = sort(round.(norm.(pos_default .- Ref(pos)), digits=11))
            radii_rotated = sort(round.(norm.(pos_rotated .- Ref(pos)), digits=11))
            @test radii_default ≈ radii_rotated
            @test !all(pos_default .≈ pos_rotated)

            # the default is reproducible, and passing the default basis reproduces it
            @test all(pos_default .≈ position.(first.(BMO.rays.(BMO.beams(
                CollimatedSource(pos, dir, diameter; num_rays, num_rings))))))
            @test all(pos_default .≈ position.(first.(BMO.rays.(BMO.beams(
                CollimatedSource(pos, dir, diameter; num_rays, num_rings, basis = b0))))))

            # a basis component along dir is projected out, so it changes nothing
            @test all(pos_default .≈ position.(first.(BMO.rays.(BMO.beams(
                CollimatedSource(pos, dir, diameter; num_rays, num_rings,
                    basis = b0 + 5 * dir))))))
        end

        @testset "Testing non-unit dir" begin
            # `rotate3d` is only a rotation for a unit-length axis, so a non-unit `dir`
            # used to scale `helper` on every step and smear the rings into a spiral
            scaled = CollimatedSource(pos, 3dir, diameter; num_rays, num_rings)
            radii = norm.(position.(first.(BMO.rays.(BMO.beams(scaled)))) .- Ref(pos))
            @test maximum(radii) ≈ diameter / 2
            @test length(unique(round.(radii, digits=11))) == num_rings
            @test BMO.direction(scaled) ≈ normalize(dir)
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException CollimatedSource(pos, dir, diameter; num_rays=100, num_rings=10)
            # basis parallel to dir has no component in the pupil plane
            @test_throws ErrorException CollimatedSource(pos, dir, diameter; num_rays, num_rings,
                basis = dir)
        end
    end

    @testset "Uniform disc source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 1]
        diameter = 2BMO.inch
        num_rays = 500

        source = UniformDiscSource(pos, dir, diameter, lambda; num_rays)

        @testset "Testing basis kwarg" begin
            dir_n = normalize(dir)
            b0 = BMO.normal3d(dir_n)
            rotated = UniformDiscSource(pos, dir, diameter, lambda; num_rays,
                basis = BMO.rotate3d(dir_n, deg2rad(30)) * b0)

            pos_default = position.(first.(BMO.rays.(BMO.beams(source))))
            pos_rotated = position.(first.(BMO.rays.(BMO.beams(rotated))))

            # rotating the basis spins the sunflower pattern about its own axis: the
            # equal-area radial distribution is invariant, the ray positions are not
            radii_default = sort(round.(norm.(pos_default .- Ref(pos)), digits=11))
            radii_rotated = sort(round.(norm.(pos_rotated .- Ref(pos)), digits=11))
            @test radii_default ≈ radii_rotated
            @test !all(pos_default .≈ pos_rotated)
            # the rotated pattern stays in the pupil plane
            @test all(isapprox.(dot.(pos_rotated .- Ref(pos), Ref(dir_n)), 0, atol = 1e-12))

            # the default is reproducible, and passing the default basis reproduces it
            @test all(pos_default .≈ position.(first.(BMO.rays.(BMO.beams(
                UniformDiscSource(pos, dir, diameter, lambda; num_rays, basis = b0))))))
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException UniformDiscSource(pos, dir, diameter, lambda;
                num_rays, basis = dir)
        end
    end

    @testset "Uniform point source" begin
        # define parameters
        lambda = 486.0e-9
        pos = [0, -0.5, 0]
        dir = [0, 1, 1]
        dir_n = normalize(dir)
        θ = deg2rad(20)
        num_rays = 500

        source = UniformPointSource(pos, dir, θ, lambda; num_rays)
        start_dirs(s) = BMO.direction.(first.(BMO.rays.(BMO.beams(s))))

        @testset "Testing getters and sampling bounds" begin
            @test source isa PointSource
            @test length(source) == num_rays
            @test BMO.numerical_aperture(source) == sin(θ)
            @test position(source) == pos
            @test BMO.direction(source) ≈ dir_n
            @test BMO.wavelength(source) == lambda
            @test all(position(b) == pos for b in BMO.beams(source))
            @test all(BMO.angle3d(dir_n, d) ≤ θ + 1e-12 for d in start_dirs(source))
            # default wavelength and number of rays
            @test length(UniformPointSource(pos, dir, θ)) == 1000
            # single ray
            @test length(UniformPointSource(pos, dir, θ; num_rays = 1)) == 1
        end

        @testset "Testing equal solid angle" begin
            N = 2000
            us = UniformPointSource(pos, dir, θ; num_rays = N)
            ϑs = BMO.angle3d.(Ref(dir_n), start_dirs(us))
            frac = count(≤(θ / 2), ϑs) / N
            @test abs(frac - (1 - cos(θ / 2)) / (1 - cos(θ))) ≤ 1 / N
        end

        @testset "Testing basis kwarg" begin
            b0 = BMO.normal3d(dir_n)
            α = deg2rad(30)
            Rα = BMO.rotate3d(dir_n, α)
            rotated = UniformPointSource(pos, dir, θ, lambda; num_rays, basis = Rα * b0)

            dirs_default = start_dirs(source)
            dirs_rotated = start_dirs(rotated)
            # rotating the basis spins the sunflower pattern about its own axis
            @test all(norm(Rα * d0 - d1) < 1e-12 for (d0, d1) in zip(dirs_default, dirs_rotated))
            @test !all(dirs_default .≈ dirs_rotated)

            # the default is reproducible, and passing the default basis reproduces it
            @test all(dirs_default .== start_dirs(UniformPointSource(pos, dir, θ, lambda; num_rays)))
            @test all(dirs_default .≈ start_dirs(UniformPointSource(pos, dir, θ, lambda; num_rays, basis = b0)))

            # orientation: local x is the projected basis, local y is dir
            b = [1.0, 2.0, 0.5]
            ob = UniformPointSource(pos, dir, θ; num_rays = 10, basis = b)
            O = orientation(ob)
            @test norm(O[:, 1] - normalize(b - dot(b, dir_n) * dir_n)) < 1e-12
            @test BMO.direction(ob) == O[:, 2]
            @test norm(O' * O - I) < 1e-12
            @test abs(det(O) - 1) < 1e-12
        end

        @testset "Testing throw errors" begin
            @test_throws ErrorException UniformPointSource(pos, dir, 1.1 * π, lambda; num_rays)
            @test_throws ErrorException UniformPointSource(pos, dir, π, lambda; num_rays)
            @test_throws ErrorException UniformPointSource(pos, dir, θ, lambda; num_rays = 0)
            # basis parallel to dir has no component in the sampling plane
            @test_throws ErrorException UniformPointSource(pos, dir, θ, lambda; num_rays, basis = dir)
            @test_throws ErrorException UniformPointSource(pos, dir, θ, lambda; num_rays, basis = [0, 0, 0])
        end
    end

    @testset "set_num_rays!" begin
        pos, dir, λ = [0.1, -0.2, 0.3], normalize([0.2, 1.0, -0.1]), 633e-9
        start(b) = Vector(BMO.position(first(rays(b))))
        heading(b) = Vector(BMO.direction(first(rays(b))))
        # A regenerated source equals a new one in its current pose
        function same(a, b)
            length(a) == length(b) || return false
            return all(zip(BMO.beams(a), BMO.beams(b))) do (x, y)
                isapprox(start(x), start(y); atol = 1e-12) && isapprox(heading(x), heading(y); atol = 1e-12)
            end
        end
        # position, direction, number of rays, sampling basis
        sources = [
            (p, d, n, b) -> CollimatedSource(p, d, 5e-3, λ; num_rings = 4, num_rays = n, basis = b),
            (p, d, n, b) -> UniformDiscSource(p, d, 5e-3, λ; num_rays = n, basis = b),
            (p, d, n, b) -> PointSource(p, d, deg2rad(100), λ; num_rings = 4, num_rays = n, basis = b),
            (p, d, n, b) -> UniformPointSource(p, d, 0.3, λ; num_rays = n, basis = b)
        ]
        for make in sources
            src = make(pos, dir, 200, nothing)
            # moved and rotated, i.e. not in the pose of the constructor
            translate3d!(src, [0.01, 0.02, -0.03])
            rotate3d!(src, normalize([1.0, 0.3, 0.2]), 0.4)
            @test set_num_rays!(src, 500) === src
            M = src.orientation
            fresh = make(Vector(position(src)), Vector(M[:, 2]), 500, Vector(M[:, 1]))
            @test same(src, fresh)
            @test length(src) == 500
            @test all(b -> BMO.wavelength(first(rays(b))) == λ, BMO.beams(src))
        end
        # a point source wider than 90°, whose NA does not tell its half angle
        ps = PointSource(pos, dir, deg2rad(100), λ; num_rings = 4, num_rays = 200)
        set_num_rays!(ps, 100)
        @test maximum(b -> acosd(clamp(dot(heading(b), dir), -1, 1)), BMO.beams(ps)) ≈ 100
        # too few rays for the rings, like the constructor
        cs = CollimatedSource(pos, dir, 5e-3, λ; num_rings = 4, num_rays = 200)
        @test_throws ErrorException set_num_rays!(cs, 79)
        @test length(cs) == 200
        # given beams can not be regenerated
        wrapped = CollimatedSource(BMO.beams(cs), 5e-3, pos, dir)
        @test_throws ArgumentError set_num_rays!(wrapped, 100)
        @test_throws ArgumentError set_num_rays!(PointSource(BMO.beams(ps), 0.5, pos, dir), 100)

        # the rays start at the emission point, not at a pivot moved by set_pivot3d!
        for make in sources
            src = make(pos, dir, 200, nothing)
            rotate3d!(src, normalize([1.0, 0.3, 0.2]), 0.4)
            M = src.orientation
            set_pivot3d!(src, pos + [0.0, -0.1, 0.02])
            set_num_rays!(src, 300)
            fresh = make(pos, Vector(M[:, 2]), 300, Vector(M[:, 1]))
            @test same(src, fresh)
            @test position(src) ≈ pos + [0.0, -0.1, 0.02]
        end

        # min_num_rays is the lower bound of set_num_rays!, the source keeps its rays below it
        for (make, n_min) in zip(sources, (80, 1, 80, 1))
            src = make(pos, dir, 200, nothing)
            @test BMO.min_num_rays(src) == n_min
            @test_throws ErrorException set_num_rays!(src, n_min - 1)
            @test length(src) == 200
            @test length(set_num_rays!(src, n_min)) == n_min
        end
        @test isnothing(BMO.min_num_rays(wrapped))
        @test isnothing(BMO.min_num_rays(Beam(pos, dir, λ)))
    end

    @testset "wavelength" begin
        λ = 633e-9
        beam = Beam([0, 0, 0], [0, 1, 0], λ)
        @test BMO.wavelength(beam) == λ
        @test BMO.wavelength(first(rays(beam))) == λ
        @test BMO.wavelength(UniformDiscSource([0, 0, 0], [0, 1, 0], 1e-3, λ; num_rays = 10)) == λ
        @test BMO.wavelength(PointSource([0, 0, 0], [0, 1, 0], 0.1, λ; num_rings = 2, num_rays = 50)) == λ
        @test Base.ispublic(BeamletOptics, :wavelength)
    end

    @testset "empty!" begin
        mirror = RoundPlanoMirror(25e-3, 5e-3)
        translate3d!(mirror, [0, 0.1, 0])
        source = CollimatedSource([0, 0, 0], [0, 1, 0], 5e-3, 1e-6; num_rings = 2)
        pos, R = position(source), BMO.orientation(source)
        solve_system!(System([mirror]), source)
        @test all(b -> length(rays(b)) == 2, BMO.beams(source))
        @test empty!(source) === source
        @test all(b -> length(rays(b)) == 1, BMO.beams(source))
        @test all(b -> isnothing(BMO.intersection(first(rays(b)))), BMO.beams(source))
        @test position(source) == pos
        @test BMO.orientation(source) == R
    end
end

end # MODULE