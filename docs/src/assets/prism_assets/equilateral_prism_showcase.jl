using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

# N-SF11 dense flint, Schott Sellmeier coefficients
NSF11 = SellmeierEquation(1.73759695, 0.313747346, 1.89878101, 0.013188707, 0.0623068142, 155.23629)

side = 25e-3
prism = EquilateralPrism(side, side, NSF11)
system = StaticSystem([prism])

# Minimum deviation at 550 nm: the ray inside runs parallel to the base
α = deg2rad(60)
θ1 = asin(NSF11(550e-9) * sin(α / 2))
ϕ = θ1 - α / 2                                  # angle of incidence above the +x-axis
r = side * sqrt(3) / 6                          # inradius
hit = [-side / 4, r / 2, 0]                     # center of the left face
dir = [cos(ϕ), sin(ϕ), 0]
pos = hit - 30e-3 * dir

λs = range(450e-9, 650e-9, length = 5)
beams = [Beam(pos, dir, λ) for λ in λs]
for beam in beams
    solve_system!(system, beam)
end

fig = Figure(size = (600, 300))
ax = Axis3(fig[1, 1], aspect = :data, azimuth = -π / 2, elevation = π / 2,
    limits = (-0.045, 0.045, -0.03, 0.015, -0.02, 0.02))
hidedecorations!(ax)
hidespines!(ax)

render!(ax, system)
for (beam, λ) in zip(beams, λs)
    render!(ax, beam; color = RGBf(wavelength_color(λ)...), flen = 40e-3)
end
# White light before the prism, lifted towards the camera above the colored rays
lift = [0, 0, 1e-3]
lines!(ax, [Point3f(pos + lift), Point3f(hit + lift)]; color = :gray40, linewidth = 2)

save("equilateral_prism_showcase.png", fig; px_per_unit = 4, update = false)
