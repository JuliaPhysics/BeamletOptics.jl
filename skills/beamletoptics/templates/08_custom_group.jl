# An own assembly type: a periscope with named mirrors that moves and traces like an ObjectGroup.
# Its mirrors do not touch, so each one interacts on its own and the type needs no interact3d.
# Parts that touch and interact as one (cemented, a Hint across a shared face) need an
# AbstractObject with the MultiShape trait and its own interact3d instead, like CubeBeamsplitter.
using BeamletOptics
using LinearAlgebra: I
const BMO = BeamletOptics

const mm = 1e-3

# A subtype of AbstractObjectGroup needs the fields `objects`, `center` (pivot) and `dir`
# (orientation). It is mutable, since moving the group changes `center` and `dir`.
mutable struct Periscope{T, L, U} <: BMO.AbstractObjectGroup{T}
    center::BMO.Point3{T}
    dir::BMO.SMatrix{3, 3, T, 9}
    const objects::Tuple{L, U}
end

# Named parts
lower_mirror(p::Periscope) = p.objects[1]
upper_mirror(p::Periscope) = p.objects[2]

# Built along +y: the lower mirror folds +y -> +z, the upper one +z -> +y, `height` above it
function Periscope(height::Real; diameter::Real = 25.4mm)
    lower = RoundPlanoMirror(diameter, 6mm)
    xrotate3d!(lower, deg2rad(-45))
    upper = RoundPlanoMirror(diameter, 6mm)
    xrotate3d!(upper, deg2rad(135))
    translate3d!(upper, [0, 0, height])
    objects = (lower, upper)
    # the members of a group are either all static or all movable
    allequal(BMO.is_static, objects) || throw(ArgumentError("mixed static and movable parts"))
    # pivot halfway up, frame of the global axes
    return Periscope(BMO.Point3(0.0, 0.0, height / 2), BMO.SMatrix{3, 3}(1.0I), objects)
end

periscope = Periscope(100mm)
translate3d!(periscope, [0, 50mm, 0])        # moves both mirrors and the pivot

det = Detector(20mm)
translate3d!(det, [0, 250mm, 100mm])

system = System([periscope, det])            # a group is flattened, each mirror is traced on its own
beam = Beam([0, 0, 0], [0, 1, 0], 633e-9)
solve_system!(system, beam)
println("upper mirror at ", position(upper_mirror(periscope)), ", ",
    length(rays(beam)), " ray segments")

# Turn the whole periscope about its pivot and trace again: the parallel mirrors keep the
# beam direction (+y), only the beam position on the detector changes
zrotate3d!(periscope, deg2rad(2))
solve_system!(system, beam; initialize = true)
println("after the turn: end direction ", direction(last(rays(beam))))
