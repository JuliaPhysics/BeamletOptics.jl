# Mirrors

Mirrors are perfect reflectors (R = 1, no losses). A `Mirror` reflects at **every** surface crossing,
including its back and edges. Curved mirrors can have a central hole via `hole_diameter`.

## Plano

| Constructor | Geometry at spawn |
|-------------|-------------------|
| `RoundPlanoMirror(diameter, thickness; hole_diameter = nothing)` | reflecting face at y = 0, normal along y, substrate towards +y |
| `SquarePlanoMirror(width, thickness)` | as above, square |
| `RectangularPlanoMirror(width_x, height_z, thickness_y)` | as above, rectangular |
| `SquarePlanoMirror2D(edge_length)` | zero-thickness square in the x-z plane |
| `RightAnglePrismMirror(leg_length, height)` | right-angle prism with a reflecting hypotenuse; legs in x/y, height in z |

To fold by 90°, rotate by 45° about the axis perpendicular to the fold plane:
`zrotate3d!(m, deg2rad(45))` sends a +y beam to +x, and `xrotate3d!(m, deg2rad(-45))` sends it to +z.

## Curved

All curved mirrors have their vertex at the origin and open towards **−y**: the concave side faces a
beam that comes from negative y and travels along +y.

| Constructor | Notes |
|-------------|-------|
| `SphericalMirror(radius, thickness, diameter; hole_diameter)` | concave sphere |
| `ParabolicMirror(f, diameter; thickness, hole_diameter)` | focus at (0, −f, 0) |
| `OffAxisParabolicMirror(rfl, diameter; angle = 90 (deg), thickness, hole_diameter, hole_axis = :collimated)` | `rfl` = reflected focal length; `hole_axis = :collimated` or `:focused` |
| `ConicMirror(R, k, diameter; thickness, hole_diameter)` | R > 0 concave; k = −1 parabola, k = 0 sphere |
| `OffAxisConicMirror(R, k, x_off, diameter; …)` | off-axis section at `x_off` |
| `EllipsoidalMirror(s, s′, diameter; …)` | conjugate distances with the same sign (positive towards −y) |
| `OffAxisEllipsoidalMirror(s, s′, x_off, diameter; …)` | |
| `HyperbolicMirror(s, s′, diameter; …)` | conjugate distances with opposite signs |
| `OffAxisHyperbolicMirror(s, s′, x_off, diameter; …)` | |

Before building a setup around an off-axis mirror, check its exact frame with
`scripts/api_lookup.jl OffAxisParabolicMirror` and a quick render.

## Retroreflector

`Retroreflector(scale)` is a corner-cube retroreflector mesh; `scale = 1e-3` gives mm dimensions.
