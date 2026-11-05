using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

# N-BK7, Schott Sellmeier coefficients
NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)

aperture = 10e-3
prism = DovePrism(4.2aperture, aperture, aperture, NBK7)
system = StaticSystem([prism])

# Three parallel rays along +y, the order across the prism is reversed after it
offsets = (-0.25aperture, 0.0, 0.25aperture)
colors = (:red, :green2, :blue)
beams = [Beam([x, -0.04, 0], [0, 1, 0], 532e-9) for x in offsets]
for beam in beams
    solve_system!(system, beam)
end

fig = Figure(size = (600, 240))
ax = Axis3(fig[1, 1], aspect = :data, azimuth = 0.0, elevation = π / 2,
    limits = (-0.01, 0.01, -0.045, 0.045, -0.01, 0.01))
hidedecorations!(ax)
hidespines!(ax)

render!(ax, system)
for (beam, color) in zip(beams, colors)
    render!(ax, beam; color = color, flen = 20e-3)
end

save("dove_prism_showcase.png", fig; px_per_unit = 4, update = false)
