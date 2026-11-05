using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

# N-BK7, Schott Sellmeier coefficients
NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)

leg = 20e-3
prism = RightAnglePrism(leg, leg, NBK7)
system = StaticSystem([prism])

# Rays along +y enter the leg with the normal -y and are totally reflected at the hypotenuse
offsets = (-0.3leg, -0.15leg, 0.0)
colors = (:red, :green2, :blue)
beams = [Beam([x, -0.03, 0], [0, 1, 0], 532e-9) for x in offsets]
for beam in beams
    solve_system!(system, beam)
end

fig = Figure(size = (600, 400))
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -π / 2, elevation = π / 2,
    limits = (-0.04, 0.02, -0.03, 0.01, -0.02, 0.02))
hidedecorations!(ax)
hidespines!(ax)

render!(ax, system)
for (beam, color) in zip(beams, colors)
    render!(ax, beam; color = color, flen = 20e-3)
end

save("right_angle_prism_showcase.png", fig; px_per_unit = 4, update = false)
