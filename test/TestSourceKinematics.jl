module TestSourceKinematics

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

"Lens followed by a 45° mirror, beams start at y = -50 mm and travel along +y."
function kinematics_test_system()
    lens = SphericalLens(50mm, -50mm, 5mm, 25mm)
    mirror = RoundPlanoMirror(25mm, 5mm)
    zrotate3d!(mirror, deg2rad(45))
    translate3d!(mirror, [0, 60mm, 0])
    return System([lens, mirror])
end

"Tests that `b` and all of its sub-beams are in the untraced start state."
function is_reset(b::Beam)
    return length(BMO.rays(b)) == 1 &&
           isnothing(BMO.intersection(first(BMO.rays(b)))) &&
           isempty(BMO.children(b))
end
is_reset(g::GaussianBeamlet) = all(is_reset, (g.chief, g.waist, g.divergence)) && isempty(BMO.children(g))
is_reset(a::AstigmaticGaussianBeamlet) = all(is_reset, BMO._component_beams(a)) && isempty(BMO.children(a))
is_reset(bg::BMO.AbstractBeamGroup) = all(is_reset, BMO.beams(bg))

end_point(b::Beam) = position(last(BMO.rays(b)))

start_rays(b) = map(BMO.first_ray, BMO._component_beams(b))

@testset "Source kinematics" begin
    system = kinematics_test_system()
    R = BMO.rotate3d(normalize([1.0, 2.0, 3.0]), 0.7)

    @testset "empty! resets traced beams" begin
        b = Beam([1mm, -50mm, 0], [0, 1, 0])
        solve_system!(system, b)
        @test length(BMO.rays(b)) > 1
        @test empty!(b) === b
        @test is_reset(b)

        g = GaussianBeamlet([0, -50mm, 0], [0, 1, 0], 1e-6, 1mm)
        solve_system!(system, g)
        @test length(BMO.rays(g.chief)) > 1
        @test empty!(g) === g
        @test is_reset(g)

        a = AstigmaticGaussianBeamlet([0, -50mm, 0], [0, 1, 0], 1e-6, 1mm)
        solve_system!(system, a)
        @test length(BMO.rays(a.c)) > 1
        @test empty!(a) === a
        @test is_reset(a)
    end

    @testset "Every verb resets a traced source" begin
        verbs = (
            s -> translate3d!(s, [1mm, 0, 0]),
            s -> translate_to3d!(s, [1mm, -50mm, 1mm]),
            s -> rotate3d!(s, BMO.rotate3d([0, 0, 1.0], 1e-3)),
            s -> rotate3d!(s, [0, 0, 1.0], 1e-3),
            s -> xrotate3d!(s, 1e-3),
            s -> yrotate3d!(s, 1e-3),
            s -> zrotate3d!(s, 1e-3),
            s -> align3d!(s, normalize([1e-3, 1, 0]))
        )
        for verb in verbs
            b = Beam([1mm, -50mm, 0], [0, 1, 0])
            solve_system!(system, b)
            @test verb(b) === nothing
            @test is_reset(b)

            g = GaussianBeamlet([0, -50mm, 0], [0, 1, 0], 1e-6, 1mm)
            solve_system!(system, g)
            verb(g)
            @test is_reset(g)

            a = AstigmaticGaussianBeamlet([0, -50mm, 0], [0, 1, 0], 1e-6, 1mm)
            solve_system!(system, a)
            verb(a)
            @test is_reset(a)

            cs = CollimatedSource([0, -50mm, 0], [0, 1, 0], 5mm; num_rings = 2, num_rays = 40)
            solve_system!(system, cs)
            @test verb(cs) === nothing
            @test is_reset(cs)
        end
    end

    @testset "Rotate traced Beam by R and back" begin
        b = Beam([1mm, -50mm, 0.5mm], [0, 1, 0])
        solve_system!(system, b)
        p0 = position(b)
        d0 = BMO.direction(b)
        e0 = end_point(b)
        rotate3d!(b, R)
        @test BMO.direction(b) ≈ R * d0
        rotate3d!(b, R')
        @test norm(position(b) - p0) < 1e-12
        @test norm(BMO.direction(b) - d0) < 1e-12
        solve_system!(system, b)
        @test norm(end_point(b) - e0) < 1e-12
    end

    @testset "Ray verbs" begin
        r = Ray([1.0, 2, 3], [0, 1, 0])
        BMO.intersection!(r, BMO.Intersection(1.0, [0.0, -1, 0]))
        @test translate3d!(r, [1, 1, 1]) === nothing
        @test position(r) == [2, 3, 4]
        @test isnothing(BMO.intersection(r))
        translate_to3d!(r, [0, 0, 0])
        @test position(r) == zeros(3)
        @test rotate3d!(r, R) === nothing
        @test BMO.direction(r) ≈ R * [0, 1, 0]
        @test position(r) == zeros(3)
        align3d!(r, [0, 0, 1])
        @test BMO.direction(r) ≈ [0, 0, 1]
        zrotate3d!(r, π / 2)
        @test BMO.direction(r) ≈ [0, 0, 1]
        xrotate3d!(r, π / 2)
        @test BMO.direction(r) ≈ [0, -1, 0]
        yrotate3d!(r, π / 2)
        @test BMO.direction(r) ≈ [0, -1, 0]
    end

    @testset "PolarizedRay E0 rotation" begin
        E0 = [1.0 + 0.5im, 0, 0.3im]
        r = PolarizedRay([0, 0, 0], [0, 1, 0], 1e-6, E0)
        rotate3d!(r, R)
        @test abs(dot(BMO.direction(r), BMO.polarization(r))) < 1e-10
        @test BMO.polarization(r) ≈ R * E0
        @test norm(BMO.polarization(r)) ≈ norm(E0)

        b = Beam([0, -50mm, 0], [0, 1, 0], 1e-6, [1.0, 0, 0])
        solve_system!(system, b)
        E_old = BMO.polarization(first(BMO.rays(b)))
        zrotate3d!(b, deg2rad(30))
        @test is_reset(b)
        E_new = BMO.polarization(first(BMO.rays(b)))
        @test abs(dot(BMO.direction(b), E_new)) < 1e-10
        @test E_new ≈ BMO.rotate3d([0, 0, 1.0], deg2rad(30)) * E_old
    end

    @testset "Rotation about an arbitrary pivot" begin
        p = [0.5, -1.0, 2.0]
        # ray: position rotated about p, direction and E0 rotated by R
        E0 = [1.0 + 0.5im, 0, 0.3im]
        pr = PolarizedRay([1.0, 2, 3], [0, 1, 0], 1e-6, E0)
        @test rotate3d!(pr, R, p) === nothing
        @test position(pr) ≈ p + R * ([1.0, 2, 3] - p)
        @test BMO.direction(pr) ≈ R * [0, 1, 0]
        @test BMO.polarization(pr) ≈ R * E0
        # beam group: center and all beams rotated about p, orientation rotated, traced beams reset
        cs = CollimatedSource([1.0, 0, 0], [0, 1, 0], 5mm; num_rings = 2, num_rays = 40)
        c0, M0 = position(cs), BMO.orientation(cs)
        starts = [position(b) for b in BMO.beams(cs)]
        solve_system!(system, cs)
        rotate3d!(cs, R, BMO.Point3(p...))
        @test position(cs) ≈ p + R * (c0 - p)
        @test BMO.orientation(cs) ≈ R * M0
        for (b, s) in zip(BMO.beams(cs), starts)
            @test position(b) ≈ p + R * (s - p)
        end
        @test is_reset(cs)
    end

    @testset "Ray rotation via axis and angle" begin
        axis = normalize([1.0, 2.0, 3.0])
        θ = 0.7
        r = Ray([1.0, 2, 3], [0, 1, 0])
        BMO.intersection!(r, BMO.Intersection(1.0, [0.0, -1, 0]))
        @test rotate3d!(r, axis, θ) === nothing
        @test BMO.direction(r) ≈ R * [0, 1, 0]
        @test position(r) == [1, 2, 3]
        @test isnothing(BMO.intersection(r))

        E0 = [1.0 + 0.5im, 0, 0.3im]
        pr = PolarizedRay([0, 0, 0], [0, 1, 0], 1e-6, E0)
        @test rotate3d!(pr, axis, θ) === nothing
        @test BMO.direction(pr) ≈ R * [0, 1, 0]
        @test abs(dot(BMO.direction(pr), BMO.polarization(pr))) < 1e-10
        @test BMO.polarization(pr) ≈ R * E0
    end

    @testset "AstigmaticGaussianBeamlet rotation" begin
        w = 1mm
        a = AstigmaticGaussianBeamlet([1.0, 0, 0], [0, 1, 0], 1e-6, w; support = [1.0, 0, 0])
        chief = BMO.polarization(a)
        dxp_dir = BMO.direction(first(BMO.rays(a.dxp)))
        Rz = BMO.rotate3d([0, 0, 1.0], π / 2)
        rotate3d!(a, Rz)
        c = first(BMO.rays(a.c))
        @test position(c) == [1, 0, 0]
        @test BMO.direction(c) ≈ [-1, 0, 0]
        @test abs(dot(BMO.direction(c), BMO.polarization(c))) < 1e-10
        @test BMO.polarization(c) ≈ Rz * chief
        wxp = first(BMO.rays(a.wxp))
        @test position(wxp) ≈ [1, w, 0]
        @test BMO.direction(wxp) ≈ [-1, 0, 0]
        dxp = first(BMO.rays(a.dxp))
        @test position(dxp) ≈ [1, 0, 0]
        @test BMO.direction(dxp) ≈ Rz * dxp_dir
        # arbitrary rotation keeps E0 orthogonal
        rotate3d!(a, R)
        @test abs(dot(BMO.direction(a), BMO.polarization(a))) < 1e-10
        @test BMO.polarization(a) ≈ R * Rz * chief
    end

    @testset "GaussianBeamlet rotated about its chief" begin
        pos = [1mm, -2mm, 3mm]
        g_ref = GaussianBeamlet(pos, [0, 1, 0], 1e-6, 0.5mm; support = [1.0, 0, 0])
        g_rot = GaussianBeamlet(pos, [0, 1, 0], 1e-6, 0.5mm; support = [1.0, 0, 0])
        rotate3d!(g_rot, R)
        @test position(g_rot) == position(g_ref)
        @test BMO.direction(g_rot) ≈ R * [0, 1, 0]
        # start rays are rigidly rotated about the chief start
        for (r_rot, r_ref) in zip(start_rays(g_rot), start_rays(g_ref))
            @test norm(R * (position(r_ref) - pos) + pos - position(r_rot)) < 1e-12
            @test norm(R * BMO.direction(r_ref) - BMO.direction(r_rot)) < 1e-12
        end
        # `gauss_parameters` itself is only conditioned to ≈1e-10 (relative) for
        # arbitrarily oriented beamlets, even without kinematics, hence the looser rtol
        for z in (0.1, 0.5, 2.0)
            p_ref = gauss_parameters(g_ref, z)
            p_rot = gauss_parameters(g_rot, z)
            @test all(isapprox.(p_rot, p_ref; rtol = 1e-8))
        end
    end

    @testset "CollimatedSource rotated about its center" begin
        c = [1mm, -50mm, 2mm]
        dir0 = [0, 1.0, 0]
        cs = CollimatedSource(c, dir0, 10mm; num_rings = 3, num_rays = 100)
        p0 = [position(b) for b in BMO.beams(cs)]
        axis = normalize([1.0, 0, 1])
        Rc = BMO.rotate3d(axis, π / 2)
        rotate3d!(cs, axis, π / 2)
        @test all(norm(Rc * (p0[i] - c) + c - position(b)) < 1e-12
        for (i, b) in enumerate(BMO.beams(cs)))
        @test all(BMO.direction(b) ≈ Rc * dir0 for b in BMO.beams(cs))
        @test position(cs) == c
        @test BMO.direction(cs) ≈ Rc * dir0
    end

    @testset "Group getters" begin
        pos = [1.0, 2, 3]
        us = UniformDiscSource(pos, [0, 1, 0], 10mm; num_rays = 50)
        @test position(us) == pos
        @test position(first(BMO.beams(us))) != pos
        @test BMO.direction(us) == [0, 1, 0]
        # positional ctor takes explicit pos/dir, dir is normalized
        cs = CollimatedSource(BMO.beams(us), 10mm, pos, [0, 2, 0])
        @test position(cs) == pos
        @test BMO.direction(cs) == [0, 1, 0]
        # setters
        BMO.position!(cs, [0, 0, 1])
        M = BMO.rotate3d([1.0, 0, 0], π / 2)
        BMO.orientation!(cs, M)
        @test position(cs) == [0, 0, 1]
        @test orientation(cs) == M
        @test BMO.direction(cs) == M[:, 2]
        @test BMO.direction(cs) ≈ [0, 0, 1]
        # positional ctor with explicit orientation, validated
        cm = CollimatedSource(BMO.beams(us), 10mm, pos, M)
        @test orientation(cm) == M
        @test_throws ArgumentError CollimatedSource(BMO.beams(us), 10mm, pos, 2M)
        @test_throws ArgumentError CollimatedSource(BMO.beams(us), 10mm, pos, -M)
        @test_throws ArgumentError CollimatedSource(BMO.beams(us), 10mm, pos, M[:, [2, 1, 3]])
        @test_throws ArgumentError PointSource(BMO.beams(us), 0.1, pos, -M)
        @test_throws ArgumentError AstigmaticBeamGroup(
            BMO.beams(BMO.CollimatedGaussianBeamletSource(pos, [0, 1, 0], 4mm, 1e-6, 1mm; n_grid = 2)),
            pos, -M)
        # beamlet sources store their pos/dir
        ag = BMO.CollimatedGaussianBeamletSource(pos, [0, 0, 2], 4mm, 1e-6, 1mm; n_grid = 4)
        @test position(ag) == pos
        @test BMO.direction(ag) == [0, 0, 1]
    end

    @testset "Group verbs" begin
        pos = [0, -50mm, 0]
        ps = PointSource(pos, [0, 1, 0], deg2rad(2); num_rings = 2, num_rays = 40)
        d0 = [BMO.direction(b) for b in BMO.beams(ps)]
        translate3d!(ps, [1mm, 0, 0])
        @test position(ps) ≈ pos + [1mm, 0, 0]
        @test all(position(b) ≈ pos + [1mm, 0, 0] for b in BMO.beams(ps))
        translate_to3d!(ps, [0, 0, 0])
        @test position(ps) ≈ zeros(3)
        @test all(norm(position(b)) < 1e-15 for b in BMO.beams(ps))
        zrotate3d!(ps, π / 2)
        Rz = BMO.rotate3d([0, 0, 1.0], π / 2)
        @test BMO.direction(ps) ≈ [-1, 0, 0]
        @test all(BMO.direction(b) ≈ Rz * d0[i] for (i, b) in enumerate(BMO.beams(ps)))
        xrotate3d!(ps, π / 2)
        yrotate3d!(ps, π / 2)
        align3d!(ps, [0, 0, 1])
        @test BMO.direction(ps) ≈ [0, 0, 1]
        @test BMO.direction(first(BMO.beams(ps))) ≈ [0, 0, 1]
        @test position(ps) ≈ zeros(3) atol = 1e-15

        # astigmatic beamlet group, pivot is the group center
        ag = BMO.CollimatedGaussianBeamletSource([0, 0, 0], [0, 1, 0], 4mm, 1e-6, 1mm; n_grid = 3)
        p0 = [position(b) for b in BMO.beams(ag)]
        zrotate3d!(ag, π / 2)
        @test all(position(b) ≈ Rz * p0[i] for (i, b) in enumerate(BMO.beams(ag)))
        @test all(abs(dot(BMO.direction(b), BMO.polarization(b))) < 1e-10 for b in BMO.beams(ag))
        @test BMO.direction(ag) ≈ [-1, 0, 0]
    end

    @testset "Group orientation of every source constructor" begin
        pos = [1mm, -2mm, 3mm]
        dir = [0.3, 1.0, -0.2]
        dn = normalize(dir)
        b = [1.0, 0.4, 2.0]
        b_proj = normalize(b - dot(b, dn) * dn)
        # second grid axis for tuple bases, only the first one sets the orientation
        b2 = cross(dn, b_proj)
        x = collect(range(-1mm, 1mm, length = 5))
        amp = ones(5, 5)
        phase = zeros(5, 5)
        is_valid(O) = norm(O' * O - I) < 1e-12 && abs(det(O) - 1) < 1e-12
        # sources without basis kwarg or with the default basis
        groups = Any[
            PointSource(pos, dir, 0.1; num_rings = 2, num_rays = 40),
            UniformPointSource(pos, dir, 0.1; num_rays = 20),
            CollimatedSource(pos, dir, 5mm; num_rings = 2, num_rays = 40),
            UniformDiscSource(pos, dir, 5mm; num_rays = 20),
            CollimatedGaussianBeamletSource(pos, dir, 4mm, 1e-6, 1mm; n_grid = 2),
            SphericalGaussianBeamletSource(pos, dir, 0.1, 1e-6; num_rings = 2, num_rays = 40),
            EllipticalGaussianBeamletSource(pos, dir, 0.1, 0.05, 1e-6; num_rings = 2, num_rays = 40),
            GaussianBeamletDecomposition(pos, dir, 1e-6, 1mm; n_grid = 3),
            WavefrontBeamletDecomposition(x, x, amp, phase, dir, 1e-6)
        ]
        # sources with basis kwarg
        based = Any[
            PointSource(pos, dir, 0.1; num_rings = 2, num_rays = 40, basis = b),
            UniformPointSource(pos, dir, 0.1; num_rays = 20, basis = b),
            CollimatedSource(pos, dir, 5mm; num_rings = 2, num_rays = 40, basis = b),
            UniformDiscSource(pos, dir, 5mm; num_rays = 20, basis = b),
            CollimatedGaussianBeamletSource(pos, dir, 4mm, 1e-6, 1mm; n_grid = 2, basis = (b, b2)),
            SphericalGaussianBeamletSource(pos, dir, 0.1, 1e-6; num_rings = 2, num_rays = 40, basis = b),
            EllipticalGaussianBeamletSource(pos, dir, 0.1, 0.05, 1e-6; num_rings = 2, num_rays = 40, basis = b),
            GaussianBeamletDecomposition(pos, dir, 1e-6, 1mm; n_grid = 3, basis = (b, b2)),
            WavefrontBeamletDecomposition(x, x, amp, phase, dir, 1e-6; basis = (b, b2))
        ]
        for bg in vcat(groups, based)
            O = orientation(bg)
            @test O isa BMO.SMatrix{3, 3, Float64, 9}
            @test is_valid(O)
            @test BMO.direction(bg) == O[:, 2]
            @test norm(O[:, 2] - dn) < 1e-12
        end
        for bg in based
            @test norm(orientation(bg)[:, 1] - b_proj) < 1e-12
        end
        # default basis is the deterministic normal of dir
        for bg in groups
            @test norm(orientation(bg)[:, 1] - BMO.normal3d(dn)) < 1e-12
        end
        # positional ctors with a direction vector
        us = first(groups)
        for bg in (PointSource(BMO.beams(us), 0.1, pos, 2dir),
            CollimatedSource(BMO.beams(us), 5mm, pos, 2dir),
            AstigmaticBeamGroup(BMO.beams(groups[5]), pos, 2dir))
            @test is_valid(orientation(bg))
            @test norm(BMO.direction(bg) - dn) < 1e-12
        end
    end

    @testset "Group roll and orientation round trip" begin
        dir = normalize([0.2, 1.0, 0.1])
        cs = CollimatedSource([1mm, -50mm, 2mm], dir, 5mm; num_rings = 2, num_rays = 40,
            basis = [1.0, 0, 0])
        O0 = orientation(cs)
        α = 0.4
        rotate3d!(cs, dir, α)
        O1 = orientation(cs)
        @test norm(O1[:, 1] - BMO.rotate3d(dir, α) * O0[:, 1]) < 1e-12
        @test norm(BMO.direction(cs) - O0[:, 2]) < 1e-12
        @test norm(O1' * O1 - I) < 1e-12
        # rotate by R and back restores the orientation
        rotate3d!(cs, R)
        @test orientation(cs) ≈ R * O1
        rotate3d!(cs, R')
        @test norm(orientation(cs) - O1) < 1e-12
        # align3d! keeps the local y-axis on the target
        align3d!(cs, [0, 0, 1.0])
        @test norm(BMO.direction(cs) - [0, 0, 1]) < 1e-12
        @test norm(orientation(cs)' * orientation(cs) - I) < 1e-12
    end

    @testset "Group reset" begin
        D = 5mm
        c0 = [1mm, -50mm, 2mm]
        dir = normalize([0.3, 1.0, -0.4])
        b = [1.0, 0.2, 0.5]
        movers = (
            s -> translate3d!(s, [2mm, -1mm, 5mm]),
            s -> rotate3d!(s, R),
            s -> zrotate3d!(s, 0.3),
            s -> translate_to3d!(s, [-3mm, 4mm, 1mm]),
            s -> rotate3d!(s, normalize([1.0, 1, 0]), 1.1),
            s -> align3d!(s, [0, 0, 1.0])
        )
        function check_same_start(bg, ref)
            @test length(bg) == length(ref)
            for (b1, b2) in zip(BMO.beams(bg), BMO.beams(ref))
                @test norm(position(b1) - position(b2)) < 1e-12
                @test norm(BMO.direction(b1) - BMO.direction(b2)) < 1e-12
            end
        end

        # CollimatedSource
        cs = CollimatedSource(c0, dir, D; num_rings = 3, num_rays = 60, basis = b)
        foreach(m -> m(cs), movers)
        solve_system!(system, cs)
        @test reset_rotation3d!(cs) === nothing
        @test orientation(cs) == I
        @test is_reset(cs)
        check_same_start(cs, CollimatedSource(position(cs), [0, 1, 0], D;
            num_rings = 3, num_rays = 60, basis = [1, 0, 0]))
        solve_system!(system, cs)
        @test reset_translation3d!(cs) === nothing
        @test position(cs) == zeros(3)
        @test is_reset(cs)
        check_same_start(cs, CollimatedSource(zeros(3), [0, 1, 0], D;
            num_rings = 3, num_rays = 60, basis = [1, 0, 0]))
        # resetting again is a no-op
        reset_rotation3d!(cs)
        reset_translation3d!(cs)
        @test orientation(cs) == I
        @test position(cs) == zeros(3)

        # AstigmaticBeamGroup, compared on the chief rays
        e1 = normalize(b - dot(b, dir) * dir)
        ag = CollimatedGaussianBeamletSource(c0, dir, 4mm, 1e-6, 1mm; n_grid = 3,
            basis = (e1, cross(dir, e1)))
        foreach(m -> m(ag), movers)
        solve_system!(system, ag)
        reset_translation3d!(ag)
        @test position(ag) == zeros(3)
        @test is_reset(ag)
        solve_system!(system, ag)
        reset_rotation3d!(ag)
        @test orientation(ag) == I
        @test is_reset(ag)
        check_same_start(ag, CollimatedGaussianBeamletSource(zeros(3), [0, 1, 0], 4mm, 1e-6, 1mm;
            n_grid = 3, basis = ([1.0, 0, 0], [0, 0, -1.0])))
    end

    @testset "set_pivot3d! (beam group)" begin
        c0 = [1mm, -50mm, 2mm]
        dir = normalize([0.2, 1.0, 0.1])
        cs = CollimatedSource(c0, dir, 5mm; num_rings = 2, num_rays = 40, basis = [1.0, 0, 0])
        solve_system!(system, cs)
        @test !is_reset(cs)

        p_before = [position(b) for b in BMO.beams(cs)]
        d_before = [BMO.direction(b) for b in BMO.beams(cs)]
        O0 = orientation(cs)

        p = [3mm, -40mm, -1mm]
        @test set_pivot3d!(cs, p) === nothing
        @test position(cs) == p
        @test orientation(cs) == O0
        # beams are untouched (exact) and the traced beam is NOT reset
        for (i, b) in enumerate(BMO.beams(cs))
            @test position(b) == p_before[i]
            @test BMO.direction(b) == d_before[i]
        end
        @test !is_reset(cs)

        # rotate3d! now rotates every beam about the new pivot p
        Rp = BMO.rotate3d(normalize([1.0, -2.0, 0.5]), 0.6)
        rotate3d!(cs, Rp)
        for (i, b) in enumerate(BMO.beams(cs))
            @test norm(position(b) - (p + Rp * (p_before[i] - p))) < 1e-12
        end
        @test is_reset(cs)

        # reset_translation3d! afterwards moves p back to the origin, members shifted by -p
        reset_translation3d!(cs)
        @test norm(position(cs)) < 1e-12
        for (i, b) in enumerate(BMO.beams(cs))
            @test norm(position(b) - Rp * (p_before[i] - p)) < 1e-12
        end
    end

    @testset "Moving a child beam throws" begin
        pbs = RectangularPlateBeamsplitter(36mm, 25mm, 1mm, n -> 1.5)
        zrotate3d!(pbs, deg2rad(45))
        bs_system = System([pbs])
        beam = Beam([0, -50mm, 0], [0, 1, 0], 1e-6)
        solve_system!(bs_system, beam; depth_max = 4)
        child = first(BMO.children(beam))
        @test BMO.isroot(beam)
        @test !BMO.isroot(child)
        @test BMO.isroot(GaussianBeamlet([0, -50mm, 0], [0, 1, 0], 1e-6, 1mm))
        @test_throws ArgumentError translate3d!(child, [1mm, 0, 0])
        @test_throws ArgumentError rotate3d!(child, R)
        @test_throws ArgumentError zrotate3d!(child, 0.1)
        # moving the root drops the children
        translate3d!(beam, [1mm, 0, 0])
        @test is_reset(beam)
    end

    @testset "Non-unit rotation axis is consistent across all kinematic types" begin
        # rotate3d!(x, axis, θ) must not depend on the length of `axis`
        axis = [1.0, -2.0, 0.5]
        θ = 0.8
        # compared state: (position, orientation/direction) for objects, start rays for sources
        state(x::Union{BMO.AbstractShape, BMO.AbstractObject}) = (position(x), BMO.orientation(x))
        state(r::BMO.AbstractRay) = (position(r), BMO.direction(r))
        state(b::BMO.Beam) = state(BMO.first_ray(b))
        state(b::BMO.AbstractBeam) = map(c -> state(BMO.first_ray(c)), BMO._component_beams(b))
        state(bg::BMO.AbstractBeamGroup) = (position(bg), BMO.orientation(bg), map(state, BMO.beams(bg)))
        flat(s::Tuple) = reduce(vcat, map(flat, s))
        flat(s::AbstractArray{<:Number}) = collect(Float64, vec(s))
        flat(s::AbstractArray) = reduce(vcat, map(flat, s))

        ref = RoundPlanoMirror(10mm, 2mm)
        translate3d!(ref, [1mm, 2mm, 3mm])
        subjects = Any[
            BMO.CylinderSDF(1mm, 2mm),                               # shape
            BMO.CylinderSDF(1mm, 2mm) + BMO.CylinderSDF(2mm, 1mm),   # composite SDF
            BMO.CubeMesh(1mm),                                       # mesh
            ref,                                                     # object (SingleShape)
            CubeBeamsplitter(5mm, λ -> 1.5),                         # object (MultiShape)
            ObjectGroup([RoundPlanoMirror(10mm, 2mm), CubeBeamsplitter(5mm, λ -> 1.5)]),
            Ray([1.0, 0, 0], [0, 1.0, 0]),
            PolarizedRay([1.0, 0, 0], [0, 1.0, 0], 1e-6, [0, 0, 1.0]),
            Beam([1.0, 0, 0], [0, 1.0, 0]),
            GaussianBeamlet([1.0, 0, 0], [0, 1.0, 0], 1e-6, 1mm),
            BMO.AstigmaticGaussianBeamlet([1.0, 0, 0], [0, 1.0, 0], 1e-6, 1mm),
            CollimatedSource([1.0, 0, 0], [0, 1.0, 0], 5mm; num_rings = 2, num_rays = 40),
        ]
        for x in subjects
            a = deepcopy(x)
            b = deepcopy(x)
            rotate3d!(a, normalize(axis), θ)
            rotate3d!(b, 3.7 * axis, θ)
            @test flat(state(a)) ≈ flat(state(b)) atol = 1e-12
            # the rotation is a proper rotation, not a scaled/sheared transform
            @test flat(state(a)) ≉ flat(state(x))
        end
        @test_throws ArgumentError rotate3d!(Ray([0, 0, 0], [0, 1, 0]), [0, 0, 0], θ)
        @test_throws ArgumentError rotate3d!(RoundPlanoMirror(10mm, 2mm), [0, 0, 0], θ)

        @testset "Rotation about an arbitrary pivot for all source types" begin
            # rotate3d!(x, R, pivot) and rotate3d!(x, axis, θ, pivot); shapes/objects: see TestObjectGroups.jl
            pivot = [0.3, -0.7, 1.1]
            Rp = BMO.rotate3d(axis, θ)
            frame(x::BMO.AbstractBeamGroup) = BMO.orientation(x)
            frame(x::Union{BMO.AbstractRay, BMO.AbstractBeam}) = BMO.direction(x)
            for x in filter(x -> x isa Union{BMO.AbstractRay, BMO.AbstractBeam, BMO.AbstractBeamGroup}, subjects)
                p0, f0 = position(x), frame(x)
                a = deepcopy(x)
                b = deepcopy(x)
                @test rotate3d!(a, Rp, pivot) === nothing
                @test rotate3d!(b, 3.7 * axis, θ, BMO.Point3(pivot...)) === nothing
                @test position(a) ≈ pivot + Rp * (p0 - pivot) atol = 1e-12
                @test frame(a) ≈ Rp * f0 atol = 1e-12
                @test flat(state(a)) ≈ flat(state(b)) atol = 1e-12
            end
        end
    end
end

end # MODULE
