# Keplerian 3x beam expander for a Gaussian beam; beam radius before/after.
using BeamletOptics

const mm = 1e-3
const nm = 1e-9

λ = 1064nm
n_glass = λ -> 1.5066                         # constant-index callable (λ in meters)

f1, f2 = 50mm, 150mm
# thin symmetric biconvex lenses: R = 2(n-1)f
L1 = SphericalLens(2 * 0.5066 * f1, -2 * 0.5066 * f1, 5mm, 25.4mm, n_glass)
L2 = SphericalLens(2 * 0.5066 * f2, -2 * 0.5066 * f2, 5mm, 25.4mm, n_glass)
translate3d!(L2, [0, f1 + f2 + thickness(L1), 0])  # approximate afocal spacing

system = StaticSystem([L1, L2])

w0 = 1mm
beam = GaussianBeamlet([0, -20mm, 0], [0, 1, 0], λ, w0; P0 = 1e-3)
solve_system!(system, beam)

# z is the distance along the (unfolded) beam from its start point
for z in (0.0, 20mm, 400mm, 800mm)
    w, R, ψ, w0_loc = gauss_parameters(beam, z)
    println("z = ", lpad(round(z / mm; digits = 1), 6), " mm   w = ",
            round(w / mm; digits = 3), " mm")
end
