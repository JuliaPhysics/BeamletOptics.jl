# Render a system, a ray fan and a Gaussian beamlet with Makie.
# Needs a Makie backend in the environment: CairoMakie (headless, Axis3 only) or GLMakie (LScene + camera helpers).
using CairoMakie, BeamletOptics

const mm = 1e-3

lens = SphericalLens(50mm, -50mm, 6mm, 25.4mm, λ -> 1.5)
mirror = RoundPlanoMirror(25.4mm, 5mm)
zrotate3d!(mirror, deg2rad(45))               # fold +y -> +x
translate3d!(mirror, [0, 80mm, 0])
system = System([lens, mirror])

# ray fan: one Beam per ray
fan = [Beam([x, -30mm, 0], [0, 1, 0], 532e-9) for x in range(-8mm, 8mm, length = 7)]
for b in fan
    solve_system!(system, b)
end

gauss = GaussianBeamlet([0, -30mm, 5mm], [0, 1, 0], 1064e-9, 2mm)
solve_system!(system, gauss)

fig = Figure(size = (800, 500))
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -π / 2, elevation = π / 2)  # top view (x-y plane)
render!(ax, system)
for b in fan
    render!(ax, b; color = :green, flen = 40mm)  # flen: drawn length of open-ended final rays
end
render!(ax, gauss; color = :red, flen = 40mm)

save(joinpath(@__DIR__, "07_render_system.png"), fig; px_per_unit = 2)
println("saved 07_render_system.png")
