# Folding a beam with plano mirrors and moving a sub-assembly as an ObjectGroup.
using BeamletOptics

const mm = 1e-3

# Mirror 1 at the origin: fold +y -> +z (rotate about x by 45°)
m1 = RoundPlanoMirror(25.4mm, 6mm)
xrotate3d!(m1, deg2rad(-45))

# Mirror 2 100 mm above: fold +z -> +y
m2 = RoundPlanoMirror(25.4mm, 6mm)
xrotate3d!(m2, deg2rad(135))
translate3d!(m2, [0, 0, 100mm])

periscope = ObjectGroup([m1, m2])             # pivot (center) starts at the origin
set_pivot3d!(periscope, [0, 0, 50mm])         # rotate about the assembly middle instead

det = Detector(20mm)
translate3d!(det, [0, 200mm, 100mm])

system = System([periscope, det])
beam = Beam([0, -50mm, 0], [0, 1, 0], 633e-9)

solve_system!(system, beam)
println("path 1: ", length(rays(beam)), " ray segments, end direction ", direction(last(rays(beam))))

# Nudge the whole assembly and resolve (retrace is automatic; reset the detector first)
translate3d!(periscope, [0, 5mm, 0])
empty!(det)
solve_system!(system, beam)
println("path 2: detector hits = ", length(spot_diagram(det)))
