# Malus's law with a polarized Beam and a rotating ideal polarization filter.
using BeamletOptics

const mm = 1e-3

pf = RoundPolarizationFilter(25.4mm)          # zero thickness, normal along y, transmits x
det = Detector(10mm)
translate3d!(det, [0, 50mm, 0])
system = System([pf, det])

E0 = [1.0, 0, 0]                              # V/m, must be ⟂ to the propagation direction
beam = Beam([0, -50mm, 0], [0, 1, 0], 633e-9, E0)   # 4-arg Beam -> PolarizedRay beam

for θ in 0:15:90
    θ > 0 && yrotate3d!(pf, deg2rad(15))      # rotate filter about the optical axis
    solve_system!(system, beam)
    # A fully blocked ray terminates at the filter: only the incident ray is left in the beam.
    passed = length(rays(beam)) > 1
    I_out = passed ? sum(abs2, last(rays(beam)).E0) : 0.0
    println("θ = ", lpad(θ, 2), "°   |E|² = ", round(I_out; digits = 4),
            "   cos²θ = ", round(cosd(θ)^2; digits = 4),
            "   axis = ", round.(transmission_axis(pf); digits = 3))
end
