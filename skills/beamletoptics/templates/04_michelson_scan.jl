# Michelson interferometer with a Gaussian beamlet: detector power vs. mirror displacement.
using BeamletOptics

const mm = 1e-3
const cm = 1e-2

λ = 632.8e-9
NBK7 = DiscreteRefractiveIndex([λ], [1.51509])  # exact-λ lookup, no interpolation!

# Cube beamsplitter centered at the origin; source travels along +y.
cbs = CubeBeamsplitter(25.4mm, NBK7)
zrotate3d!(cbs, deg2rad(-90))

# Arm 1 (transmitted, +y) and arm 2 (reflected, +x)
m1 = RoundPlanoMirror(25.4mm, 5mm)
translate3d!(m1, [0, 10cm, 0])
m2 = RoundPlanoMirror(25.4mm, 5mm)
zrotate3d!(m2, deg2rad(90))                   # mirror normal along x
translate3d!(m2, [10cm, 0, 0])

pd = Detector(5mm)
zrotate3d!(pd, deg2rad(90))
translate3d!(pd, [-5cm, 0, 0])                # output port on -x

system = System([cbs, m1, m2, pd])
beam = GaussianBeamlet([0, -10cm, 0], [0, 1, 0], λ, 0.5mm)

Δ = range(0, λ, length = 9)                   # move m2 by one wavelength (two fringes)
P = similar(collect(Δ))
for i in eachindex(Δ)
    i > 1 && translate3d!(m2, [step(Δ), 0, 0])
    empty!(pd)                                # detectors accumulate hits -> reset every run
    solve_system!(system, beam)
    P[i] = optical_power(pd)
end

for (d, p) in zip(Δ, P)
    println("Δx = ", lpad(round(d / 1e-9; digits = 1), 6), " nm   P = ",
            round(p / 1e-3; digits = 4), " mW")
end
