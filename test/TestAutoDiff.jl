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
end

end # MODULE
