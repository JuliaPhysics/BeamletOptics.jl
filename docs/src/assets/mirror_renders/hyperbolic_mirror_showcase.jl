using GLMakie, BeamletOptics

GLMakie.activate!(; ssao = true)

const BMO = BeamletOptics
const mm = 1e-3

## Classical Cassegrain telescope: parabolic primary with a central bore, convex hyperbolic secondary
f1 = 200mm      # primary focal length
D1 = 100mm      # primary diameter
d = 150mm       # primary-secondary vertex separation
b = 50mm        # final focus position behind the primary vertex
D2 = 30mm       # secondary diameter

# primary: vertex at the origin, opens towards -y, prime focus at (0, -f1, 0)
primary = ParabolicMirror(f1, D1; hole_diameter = 10mm, thickness=5mm)

# secondary: the prime focus lies behind it (virtual, s < 0), the final focus behind the primary (real, s′ > 0)
s = d - f1
s′ = d + b
secondary = HyperbolicMirror(s, s′, D2; thickness = 5mm)
# face the primary and move the vertex to (0, -d, 0)
zrotate3d!(secondary, π)
translate3d!(secondary, [0, -d, 0])

system = StaticSystem([primary, secondary])

## Collimated fan in the x-y plane, outside of the central obstruction by the secondary
beams = Beam{Float64, Ray{Float64}}[]
for z in [-42, -34, -26, -18, 18, 26, 34, 42] .* mm
    local beam = Beam(Ray([0, -1.25d, z], [0, 1, 0]))
    solve_system!(system, beam)
    push!(beams, beam)
end

##
cview = [
  0.315127   0.94905    5.55112e-17   0.0674507
 -0.129151   0.0428839  0.990697      0.00715206
  0.940221  -0.312195   0.136085     -0.196578
  0.0        0.0        0.0           1.0
]

fig = Figure(size = (600, 400))
display(fig)
ax = LScene(fig[1, 1], show_axis = false)

render!(ax, primary; transparency = true, alpha = 0.05)
render!(ax, secondary; transparency = true, alpha = 0.4)
for beam in beams
    # the rays leaving through the bore are untruncated, let them run a bit past the final focus
    render!(ax, beam; color = RGBAf(1, 0, 0, 0.6), show_pos = true, flen = 200mm)
end
scatter!(ax, [Point3f(0, b, 0)]; color = :black, markersize = 10)

set_view(ax, cview)

save("hyperbolic_mirror_showcase.png", fig; px_per_unit = 4, update = false)
