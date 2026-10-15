module TestAnalyticReferences

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

#=
Validation of traced results against references that do not use BMO: an exact ray trace through
spherical surfaces, closed-form results of geometrical optics, Fresnel's equations, a datasheet
and the ABCD law of Gaussian beams.
=#

"""
    exact_trace(surfaces, x0, y0; n0 = 1)

Exact meridional trace of a ray that starts at `(x0, y0)` along `+y`, independent of BMO. A surface is
`(y_vertex, R, n)`: the vertex on the `y` axis, the radius of curvature (positive if the center lies
at larger `y`, `Inf` for a plane) and the refractive index behind it. Returns the position and the
direction of the ray behind the last surface.
"""
function exact_trace(surfaces, x0, y0; n0 = 1.0)
    p = [x0, y0]
    d = [0.0, 1.0]
    n1 = n0
    for (yv, R, n2) in surfaces
        if isinf(R)
            p = p + (yv - p[2]) / d[2] * d
            N = [0.0, -1.0]
        else
            c = [0.0, yv + R]
            oc = p - c
            b = dot(oc, d)
            disc = b^2 - (dot(oc, oc) - R^2)
            # the intersection near the vertex
            p = p + (R > 0 ? -b - sqrt(disc) : -b + sqrt(disc)) * d
            N = (p - c) / R
            dot(N, d) > 0 && (N = -N)
        end
        # Snell's law in vector form, N against the ray
        μ = n1 / n2
        cosi = -dot(N, d)
        d = μ * d + (μ * cosi - sqrt(1 - μ^2 * (1 - cosi^2))) * N
        n1 = n2
    end
    return p, d
end

# height and slope of the ray (p, d) at the plane y = Y
at_plane(p, d, Y) = (p[1] + (Y - p[2]) / d[2] * d[1], d[1] / d[2])

"""
Largest difference of height [m] and slope [rad] at the plane `y = Y` between BMO and the exact
trace, for rays at the `heights` along `+y` from `y0`.
"""
function trace_errors(system, surfaces, heights; y0 = -50mm, Y = 200mm)
    err_h, err_s = 0.0, 0.0
    for h in heights
        beam = Beam([h, y0, 0.0], [0.0, 1, 0], 587.6e-9)
        solve_system!(system, beam)
        r = last(BMO.rays(beam))
        hb, sb = at_plane(BMO.position(r)[1:2], BMO.direction(r)[1:2], Y)
        he, se = at_plane(exact_trace(surfaces, h, y0)..., Y)
        err_h = max(err_h, abs(hb - he))
        err_s = max(err_s, abs(sb - se))
    end
    return err_h, err_s
end

@testset "Analytic references" begin
    @testset "Exact ray trace through spherical surfaces" begin
        n = 1.5
        hs = [0.01mm, 1mm, 4mm, 8mm, 11mm]
        lenses = [
            "bi-convex" => (SphericalLens(50mm, -50mm, 5mm, 25mm, _ -> n), [(0.0, 50mm, n), (5mm, -50mm, 1.0)]),
            "plano-convex" => (SphericalLens(50mm, Inf, 5mm, 25mm, _ -> n), [(0.0, 50mm, n), (5mm, Inf, 1.0)]),
            "convex-plano" => (SphericalLens(Inf, -50mm, 5mm, 25mm, _ -> n), [(0.0, Inf, n), (5mm, -50mm, 1.0)]),
            "bi-concave" => (SphericalLens(-50mm, 50mm, 3mm, 25mm, _ -> n), [(0.0, -50mm, n), (3mm, 50mm, 1.0)]),
            "meniscus" => (SphericalLens(30mm, 60mm, 6mm, 25mm, _ -> n), [(0.0, 30mm, n), (6mm, 60mm, 1.0)]),
            "cemented doublet" => (SphericalDoubletLens(40mm, -30mm, -120mm, 6mm, 2.5mm, 25mm, _ -> 1.62, _ -> 1.72),
                [(0.0, 40mm, 1.62), (6mm, -30mm, 1.72), (8.5mm, -120mm, 1.0)]),
        ]
        @testset "$label" for (label, (lens, surfaces)) in lenses
            err_h, err_s = trace_errors(System([lens]), surfaces, hs)
            @test err_h < 1e-9
            @test err_s < 1e-8
        end
        @testset "Sonnar 50 mm F1.5, 10 surfaces" begin
            # as in TestSonnarLens.jl
            l1 = SphericalLens(69.21mm, 433.84mm, 9.33mm, 70mm, λ -> 1.671)
            l234 = TripletLens(SphericalSurface(35.86mm, 60mm), SphericalSurface(85.87mm, 60mm),
                SphericalSurface(-646.31mm, 60mm), SphericalSurface(23.51mm, 40mm),
                11.81mm, 7.05mm, 1.9mm, λ -> 1.671, λ -> 1.4892, λ -> 1.7394)
            l567 = SphericalTripletLens(Inf, 51.09mm, -22.12mm, -103.13mm, 2.48mm, 19.81mm, 4.57mm, 42mm,
                λ -> 1.5232, λ -> 1.6578, λ -> 1.5894)
            y2 = thickness(l1) + 0.38mm
            y3 = y2 + thickness(l234) + 13.0mm + 2.24mm
            translate3d!(l234, [0, y2, 0])
            translate3d!(l567, [0, y3, 0])
            surfaces = [(0.0, 69.21mm, 1.671), (9.33mm, 433.84mm, 1.0),
                (y2, 35.86mm, 1.671), (y2 + 11.81mm, 85.87mm, 1.4892),
                (y2 + 18.86mm, -646.31mm, 1.7394), (y2 + 20.76mm, 23.51mm, 1.0),
                (y3, Inf, 1.5232), (y3 + 2.48mm, 51.09mm, 1.6578),
                (y3 + 22.29mm, -22.12mm, 1.5894), (y3 + 26.86mm, -103.13mm, 1.0)]
            err_h, err_s = trace_errors(System([l1, l234, l567]), surfaces, [0.01mm, 5mm, 15mm, 30mm]; Y = 250mm)
            @test err_h < 1e-9
            @test err_s < 1e-8
        end
    end

    @testset "Mirrors" begin
        # axis crossing of a ray at the height h, reflected by a mirror whose vertex is at y = 100 mm
        function focus_distance(mirror, h)
            translate3d!(mirror, [0, 100mm, 0])
            beam = Beam([h, 0, 0], [0.0, 1, 0], 1e-6)
            solve_system!(System([mirror]), beam)
            r = last(BMO.rays(beam))
            p, d = BMO.position(r), BMO.direction(r)
            return 100mm - (p[2] - p[1] / d[1] * d[2])
        end
        # paraxial focus of a spherical mirror at R/2
        @test focus_distance(SphericalMirror(200mm, 5mm, 50mm), 0.01mm) ≈ 100mm atol = 1e-9
        # a paraboloid focuses all heights at f
        @test all(h -> isapprox(focus_distance(ParabolicMirror(100mm, 50mm), h), 100mm; atol = 1e-10),
            (0.5mm, 5mm, 12mm, 20mm))
    end

    @testset "Minimum deviation of a prism" begin
        # n sin(A/2) = sin((A + D_min)/2) for the apex angle A
        n, A = 1.5, π / 3
        function deviation(θ)
            prism = EquilateralPrism(20mm, 10mm, n)
            zrotate3d!(prism, θ)
            translate3d!(prism, [50mm, 0, 0])
            beam = Beam([0.0, 2mm, 0], [1.0, 0, 0], 1e-6)
            solve_system!(System([prism]), beam)
            length(BMO.rays(beam)) == 3 || return NaN
            return acos(clamp(BMO.direction(last(BMO.rays(beam)))[1], -1, 1))
        end
        D_min = minimum(filter(!isnan, deviation.(range(-0.6, 0.6, 2401))))
        @test D_min ≈ 2asin(n * sin(A / 2)) - A atol = 1e-6
    end

    @testset "Sellmeier equation of N-BK7 vs datasheet" begin
        NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945,
            0.00600069867, 0.0200179144, 103.560653)
        # n_F, n_d and n_C of the Schott datasheet
        for (λ, n) in ((486.13e-9, 1.52238), (587.56e-9, 1.51680), (656.27e-9, 1.51432))
            @test NBK7(λ) ≈ n atol = 1e-5
        end
    end

    @testset "Fresnel transmission of a window" begin
        # the field behind both faces is |t12 t21| = 1 - |r|² of the incident field (Stokes relations)
        n = 1.5
        function r(θ, pol)
            c, ct = cos(θ), sqrt(1 - (sin(θ) / n)^2)
            return pol === :s ? (c - n * ct) / (c + n * ct) : (n * c - ct) / (n * c + ct)
        end
        for θ in deg2rad.((0.0, 20.0, 45.0, 70.0)), pol in (:s, :p)
            window = Lens(RectangularFlatSurface(30mm), RectangularFlatSurface(30mm), 5mm, _ -> n)
            zrotate3d!(window, θ)
            translate3d!(window, [0, 50mm, 0])
            pd = Detector(40mm)
            translate3d!(pd, [0, 120mm, 0])
            # the plane of incidence is the x-y plane
            E0 = pol === :s ? ComplexF64[0, 0, 1] : ComplexF64[1, 0, 0]
            solve_system!(System([window, pd]), Beam([0.0, 0, 0], [0.0, 1, 0], 1e-6, E0))
            @test norm(BMO.polarization(only(BMO.hits(pd)))) ≈ 1 - abs2(r(θ, pol)) atol = 1e-12
        end
    end

    @testset "Gaussian beams vs ABCD" begin
        P(d) = [1 d; 0 1]
        S(R, n1, n2) = [1 0; (n1 - n2) / (n2 * R) n1 / n2]
        prop(M, q) = (M[1, 1] * q + M[1, 2]) / (M[2, 1] * q + M[2, 2])
        radius(q, λ) = sqrt(-λ / (π * imag(1 / q)))
        λ, n, t, y_lens = 1e-6, 1.5, 5mm, 100mm
        @testset "GaussianBeamlet behind a thick lens" begin
            # a narrow beam, such that the spherical aberration of the parabasal rays is negligible
            w0 = 0.15mm
            lens = SphericalLens(50mm, -50mm, t, 25mm, _ -> n)
            translate3d!(lens, [0, y_lens, 0])
            g = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0)
            solve_system!(System([lens]), g)
            M = S(-50mm, n, 1) * P(t) * S(50mm, 1, n) * P(y_lens)
            q0 = im * π * w0^2 / λ
            for z in (20mm, 60mm, 100mm)
                w = first(BMO.gauss_parameters(g, y_lens + t + z))
                @test w ≈ radius(prop(P(z) * M, q0), λ) rtol = 1e-4
            end
        end
        @testset "AstigmaticGaussianBeamlet behind a cylindrical lens" begin
            # the curved axis focuses, the other one propagates as in a plate
            w0 = 0.5mm
            lens = Lens(CylindricalSurface(30mm, 20mm, 25mm), t, _ -> n)
            translate3d!(lens, [0, y_lens, 0])
            g = AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0)
            solve_system!(System([lens]), g)
            M_curved = S(Inf, n, 1) * P(t) * S(30mm, 1, n) * P(y_lens)
            M_flat = S(Inf, n, 1) * P(t) * S(Inf, 1, n) * P(y_lens)
            q0 = im * π * w0^2 / λ
            for z in (20mm, 50mm, 80mm)
                w1, w2 = BMO.gauss_parameters(g, y_lens + t + z)[1:2]
                w_ref = minmax(radius(prop(P(z) * M_curved, q0), λ), radius(prop(P(z) * M_flat, q0), λ))
                @test all(isapprox.(minmax(w1, w2), w_ref; rtol = 1e-3))
            end
        end
    end

    @testset "Michelson interferometer: fringe period λ/2" begin
        λ = 632.8e-9
        function power(dx)
            cm = 1e-2
            cbs = CubeBeamsplitter(25.4mm, _ -> 1.51509)
            zrotate3d!(cbs, -π / 2)
            translate3d!(cbs, [0, 23.5cm, 0])
            m1 = RoundPlanoMirror(25.4mm, 5mm)
            zrotate3d!(m1, -π / 2)
            translate3d!(m1, [5cm + dx, 23.5cm, 0])
            m2 = RoundPlanoMirror(25.4mm, 5mm)
            translate3d!(m2, [0, 30cm, 0])
            pd = Detector(8mm)
            zrotate3d!(pd, π / 2)
            translate3d!(pd, [-5cm, 23.5cm, 0])
            solve_system!(System([cbs, m1, m2, pd]), GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, 0.5mm))
            return optical_power(pd; progress = false, n = 60)
        end
        # one fringe over λ/2 of the mirror position, sampled finely enough for its contrast
        xs = range(0, λ / 2, 33)
        Ps = power.(xs)
        @test Ps[end] ≈ Ps[1] rtol = 1e-5
        @test power(λ / 2 + xs[9]) ≈ Ps[9] rtol = 1e-5
        # nearly full contrast for equal beams
        @test (maximum(Ps) - minimum(Ps)) / (maximum(Ps) + minimum(Ps)) > 0.99
    end
end

end # MODULE
