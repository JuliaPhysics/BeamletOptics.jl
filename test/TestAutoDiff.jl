module TestAutoDiff

using BeamletOptics
using ForwardDiff
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

#=
Derivatives of traced results w.r.t. the pose of a component, via ForwardDiff, compared with
central finite differences. The number type of the setup is that of the parameter `dx`, i.e. a
`ForwardDiff.Dual` within `ForwardDiff.derivative`.
=#

"""
The `x` coordinate of the hit of a ray at the height `h` on a detector, with the `component` shifted
laterally by `dx`. `component(T)` builds it in the number type `T` with its axis along `y`; a mirror
reflects the ray back towards `-y` (`mirror = true`), where the detector lies.
"""
function spot_x(component, dx; h, mirror = false)
    T = typeof(dx)
    c = component(T)
    translate3d!(c, [dx, zero(T), zero(T)])
    pd = Detector(T(40mm))
    if mirror
        zrotate3d!(pd, T(π))
        translate3d!(pd, T[0, -30mm, 0])
    else
        translate3d!(pd, T[0, 100mm, 0])
    end
    # a mirror is hit from -y, the detector lies between source and mirror and does not stop it
    y0 = mirror ? -20mm : -50mm
    beam = Beam(T[h, y0, 0], T[0, 1, 0], T(1e-6))
    system = System([c, pd])
    solve_system!(system, beam)
    return BMO.hit_point(last(BMO.hits(pd)))[1]
end

function check_derivative(component; h, mirror = false)
    f(dx) = spot_x(component, dx; h, mirror)
    δ = 1e-6
    return ForwardDiff.derivative(f, 0.0), (f(δ) - f(-δ)) / 2δ
end

n(T) = _ -> T(1.5)
lenses = [
    "bi-convex SphericalLens" => T -> SphericalLens(T(50mm), T(-50mm), T(5mm), T(25mm), n(T)),
    "bi-concave SphericalLens" => T -> SphericalLens(T(-50mm), T(50mm), T(3mm), T(25mm), n(T)),
    "meniscus SphericalLens" => T -> SphericalLens(T(30mm), T(50mm), T(4mm), T(25mm), n(T)),
    "thin meniscus (MeniscusLensSDF)" => T -> SphericalLens(T(50mm), T(30mm), T(3mm), T(25.4mm), n(T)),
    "plano-convex aspheric Lens" => T -> Lens(EvenAsphericalSurface(T(40mm), T(25mm), T(-1), T[0, 1e3]),
        CircularFlatSurface(T(25mm)), T(6mm), n(T)),
    "plano-concave aspheric Lens" => T -> Lens(EvenAsphericalSurface(T(-40mm), T(25mm), T(-1), T[0, -1e3]),
        CircularFlatSurface(T(25mm)), T(3mm), n(T)),
    "cylindrical Lens" => T -> Lens(CylindricalSurface(T(30mm), T(20mm), T(25mm)), T(5mm), n(T)),
]
mirrors = [
    "SphericalMirror" => T -> SphericalMirror(T(200mm), T(5mm), T(25mm)),
    "ParabolicMirror" => T -> ParabolicMirror(T(100mm), T(25mm)),
    "ConicMirror (hyperboloid)" => T -> ConicMirror(T(200mm), T(-2), T(25mm)),
]

pd_at(T, y) = (d = Detector(T(40mm)); translate3d!(d, T[0, y, 0]); d)
hit_x(pd) = BMO.hit_point(only(BMO.hits(pd)))[1]
bi_convex(T) = SphericalLens(T(50mm), T(-50mm), T(5mm), T(25mm), _ -> T(1.5))
# paraxial ABCD matrices in reduced angles: a translation by the reduced distance `d` (length / n),
# and a refraction or reflection with the power `P`
translation(d) = [1 d; 0 1]
refraction(P) = [1 0; -P 1]
polarized(T, x) = Beam(PolarizedRay(T[x, -50mm, 0], T[0, 1, 0], T(1e-6), Complex{T}[1, 0, 0]))

"""
A Michelson interferometer with a cube beamsplitter, the optical power on its detector for the
arm mirror `m1` shifted by `dx` along its axis.
"""
function michelson_power(dx)
    T = typeof(dx)
    cm = T(1e-2)
    cbs = CubeBeamsplitter(T(25.4mm), _ -> T(1.51509))
    zrotate3d!(cbs, T(-π / 2))
    translate3d!(cbs, T[0, 23.5cm, 0])
    m1 = RoundPlanoMirror(T(25.4mm), T(5mm))
    zrotate3d!(m1, T(-π / 2))
    translate3d!(m1, T[5cm + dx, 23.5cm, 0])
    m2 = RoundPlanoMirror(T(25.4mm), T(5mm))
    translate3d!(m2, T[0, 30cm, 0])
    pd = Detector(T(8mm))
    zrotate3d!(pd, T(π / 2))
    translate3d!(pd, T[-5cm, 23.5cm, 0])
    g = GaussianBeamlet(T[0, 0, 0], T[0, 1, 0], T(632.8e-9), T(0.5mm))
    solve_system!(System([cbs, m1, m2, pd]), g)
    return optical_power(pd, RegularGrid(num_points = 40); progress = false)
end

# label, function of the parameter, value of the parameter
parametric_cases = [
    ("radius of curvature", R -> (T = typeof(R); l = SphericalLens(R, T(-50mm), T(5mm), T(25mm), _ -> T(1.5));
        d = pd_at(T, 100mm); solve_system!(System([l, d]), Beam(T[2mm, -50mm, 0], T[0, 1, 0], T(1e-6))); hit_x(d)), 50mm),
    ("center thickness", t -> (T = typeof(t); l = SphericalLens(T(50mm), T(-50mm), t, T(25mm), _ -> T(1.5));
        d = pd_at(T, 100mm); solve_system!(System([l, d]), Beam(T[2mm, -50mm, 0], T[0, 1, 0], T(1e-6))); hit_x(d)), 5mm),
    ("lens tilt", θ -> (T = typeof(θ); l = bi_convex(T); zrotate3d!(l, θ);
        d = pd_at(T, 100mm); solve_system!(System([l, d]), Beam(T[0, -50mm, 0], T[0, 1, 0], T(1e-6))); hit_x(d)), 0.0),
    ("PolarizedRay on the axis", dx -> (T = typeof(dx); l = bi_convex(T); translate3d!(l, [dx, zero(T), zero(T)]);
        d = pd_at(T, 100mm); solve_system!(System([l, d]), polarized(T, 0)); hit_x(d)), 0.0),
    ("PolarizedRay transmitted field (Fresnel)", θ -> (T = typeof(θ); l = bi_convex(T); zrotate3d!(l, θ);
        d = pd_at(T, 100mm); solve_system!(System([l, d]), polarized(T, 0)); norm(BMO.polarization(only(BMO.hits(d))))), 0.3),
    ("AstigmaticGaussianBeamlet on the axis", dx -> (T = typeof(dx); l = bi_convex(T); translate3d!(l, [dx, zero(T), zero(T)]);
        d = pd_at(T, 100mm); g = AstigmaticGaussianBeamlet(T[0, -50mm, 0], T[0, 1, 0], T(1e-6), T(1mm));
        solve_system!(System([l, d]), g); hit_x(d)), 0.0),
    ("ThinBeamsplitter tilt, reflected beam", θ -> (T = typeof(θ); bs = ThinBeamsplitter(T(30mm)); zrotate3d!(bs, T(π / 4) + θ);
        translate3d!(bs, T[0, 50mm, 0]); d = Detector(T(40mm)); zrotate3d!(d, T(π / 2)); translate3d!(d, T[100mm, 50mm, 0]);
        solve_system!(System([bs, d]), Beam(T[0, 0, 0], T[0, 1, 0], T(1e-6))); BMO.hit_point(only(BMO.hits(d)))[2]), 0.0),
    ("CubeBeamsplitter tilt", θ -> (T = typeof(θ); c = CubeBeamsplitter(T(25mm), _ -> T(1.5)); zrotate3d!(c, θ);
        translate3d!(c, T[0, 50mm, 0]); d = pd_at(T, 150mm);
        solve_system!(System([c, d]), Beam(T[1mm, 0, 0], T[0, 1, 0], T(1e-6))); hit_x(d)), 0.05),
    ("RoundPlateBeamsplitter tilt", θ -> (T = typeof(θ); p = RoundPlateBeamsplitter(T(25mm), T(5mm), _ -> T(1.5));
        zrotate3d!(p, T(π / 4) + θ); translate3d!(p, T[0, 50mm, 0]); d = pd_at(T, 150mm);
        solve_system!(System([p, d]), Beam(T[0, 0, 0], T[0, 1, 0], T(1e-6))); hit_x(d)), 0.0),
    ("total internal reflection in a RightAnglePrism", θ -> (T = typeof(θ); p = RightAnglePrism(T(25mm), T(25mm), _ -> T(1.5));
        zrotate3d!(p, θ); translate3d!(p, T[0, 50mm, 0]); b = Beam(T[4mm, 0, 0], T[0, 1, 0], T(1e-6));
        solve_system!(System([p]), b); BMO.direction(last(BMO.rays(b)))[2]), 0.0),
    ("RoundPolarizationFilter angle (Malus)", θ -> (T = typeof(θ); f = RoundPolarizationFilter(T(25mm)); yrotate3d!(f, θ);
        translate3d!(f, T[0, 50mm, 0]); d = pd_at(T, 100mm); b = Beam(PolarizedRay(T[0, 0, 0], T[0, 1, 0], T(1e-6), Complex{T}[1, 0, 0]));
        solve_system!(System([f, d]), b); sum(abs2, BMO.polarization(only(BMO.hits(d))))), 0.5),
    ("optical_power of a clipped GaussianBeamlet", dx -> (T = typeof(dx); d = Detector(T(2mm)); translate3d!(d, [dx, T(100mm), zero(T)]);
        g = GaussianBeamlet(T[0, 0, 0], T[0, 1, 0], T(1e-6), T(1mm)); solve_system!(System([d]), g);
        optical_power(d, RegularGrid(num_points = 60, x_min = -1mm, x_max = 1mm, z_min = -1mm, z_max = 1mm); progress = false)), 0.5mm),
    ("Michelson interferometer, power vs. mirror position", michelson_power, 1e-8),
]

@testset "Automatic differentiation" begin
    @testset "$label, h = $h" for (label, c) in lenses, h in (0.0, 2mm)
        ad, fd = check_derivative(c; h)
        @test isfinite(ad)
        @test ad ≈ fd rtol = 1e-5
    end
    @testset "$label, h = $h" for (label, c) in mirrors, h in (0.0, 2mm)
        ad, fd = check_derivative(c; h, mirror = true)
        @test isfinite(ad)
        @test ad ≈ fd rtol = 1e-5
    end
    @testset "paraxial reference, $label" for (label, c, M, mirror) in [
            # front vertex at y = 0, back vertex at y = t, detector at y = 100 mm
            ("bi-convex SphericalLens", bi_convex,
                translation(100mm - 5mm) * refraction((1 - 1.5) / -50mm) * translation(5mm / 1.5) *
                refraction((1.5 - 1) / 50mm), false),
            # vertex at y = 0, detector at y = -30 mm
            ("SphericalMirror", last(mirrors[1]), translation(30mm) * refraction(2 / 200mm), true),
        ]
        # A ray on the axis and the component shifted laterally by `dx`: the ray meets it at the
        # height -dx, i.e. the spot lies at dx + A (-dx), with A the position element of the paraxial
        # ABCD matrix `M` from the vertex to the detector, exact for the derivative at dx = 0.
        ad = ForwardDiff.derivative(dx -> spot_x(c, dx; h = 0.0, mirror), 0.0)
        @test ad ≈ 1 - M[1, 1] rtol = 1e-6
    end
    @testset "refractive index" begin
        # the spot of a ray at h = 2 mm behind a bi-convex lens with the refractive index `n0`
        function spot(n0)
            T = typeof(n0)
            lens = SphericalLens(T(50mm), T(-50mm), T(5mm), T(25mm), _ -> n0)
            pd = Detector(T(40mm))
            translate3d!(pd, T[0, 100mm, 0])
            beam = Beam(T[2mm, -50mm, 0], T[0, 1, 0], T(1e-6))
            solve_system!(System([lens, pd]), beam)
            return BMO.hit_point(only(BMO.hits(pd)))[1]
        end
        ad = ForwardDiff.derivative(spot, 1.5)
        fd = (spot(1.5 + 1e-6) - spot(1.5 - 1e-6)) / 2e-6
        @test ad ≈ fd rtol = 1e-5
    end
    @testset "$label" for (label, f, x0) in parametric_cases
        ad = ForwardDiff.derivative(f, x0)
        δ = 1e-9
        fd = (f(x0 + δ) - f(x0 - δ)) / 2δ
        @test isfinite(ad)
        @test ad ≈ fd rtol = 1e-4 atol = 1e-9
    end
    @testset "GaussianBeamlet" begin
        # the hit point of a beamlet on the axis behind a lens that is shifted laterally
        function spot(dx)
            T = typeof(dx)
            lens = SphericalLens(T(50mm), T(-50mm), T(5mm), T(25mm), _ -> T(1.5))
            translate3d!(lens, [dx, zero(T), zero(T)])
            pd = Detector(T(40mm))
            translate3d!(pd, T[0, 100mm, 0])
            g = GaussianBeamlet(T[0, -50mm, 0], T[0, 1, 0], T(1e-6), T(1mm))
            solve_system!(System([lens, pd]), g)
            return BMO.hit_point(only(BMO.hits(pd)))[1]
        end
        ad = ForwardDiff.derivative(spot, 0.0)
        fd = (spot(1e-6) - spot(-1e-6)) / 2e-6
        @test ad ≈ fd rtol = 1e-5
    end
    @testset "AstigmaticGaussianBeamlet, waist offset" begin
        # the beam radius at the start with the waist offset as the only Dual, w = w0 √(1 + (z0 / zR)²)
        λ, w0, z0 = 1.064e-6, 0.5e-3, 0.2
        zR = π * w0^2 / λ
        dw = w0 * z0 / zR^2 / sqrt(1 + (z0 / zR)^2)
        pos, dir = [0.0, 0, 0], [0.0, 1, 0]
        radii(beam) = gauss_parameters(beam, 0.0)
        @test ForwardDiff.derivative(z -> radii(AstigmaticGaussianBeamlet(pos, dir, λ, w0; z0 = z))[1], z0) ≈ dw rtol = 1e-6
        @test ForwardDiff.derivative(z -> radii(AstigmaticGaussianBeamlet(pos, dir, λ, w0, w0; z0_x = z))[1], z0) ≈ dw rtol = 1e-6
        @test ForwardDiff.derivative(z -> radii(AstigmaticGaussianBeamlet(pos, dir, λ, w0, w0; z0_y = z))[2], z0) ≈ dw rtol = 1e-6
    end
    @testset "AstigmaticGaussianBeamlet, support" begin
        # the support vector [cos t, 0, sin t] as the only Dual: x of the field and of the waist ray
        λ, w0, t = 1e-6, 1e-3, 0.3
        pos, dir = [0.0, 0, 0], [0.0, 1, 0]
        beam(t) = AstigmaticGaussianBeamlet(pos, dir, λ, w0; support = [cos(t), 0, sin(t)])
        E = norm(BMO.polarization(first(BMO.rays(beam(t).c))))
        @test ForwardDiff.derivative(t -> real(BMO.polarization(first(BMO.rays(beam(t).c)))[1]), t) ≈ -E * sin(t) rtol = 1e-9
        @test ForwardDiff.derivative(t -> BMO.position(first(BMO.rays(beam(t).wxp)))[1], t) ≈ -w0 * sin(t) rtol = 1e-9
    end
end

end # MODULE
