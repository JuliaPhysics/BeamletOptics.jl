using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

##
# define spherical lenses
l1 = SphericalLens(69.21e-3, 433.84e-3, 9.33e-3, 70e-3, λ -> 1.671)
# front triplet: last surface only 40 mm clear aperture -> assembled from individual lenses
l2 = SphericalLens(35.86e-3, 85.87e-3, 11.81e-3, 60e-3, λ -> 1.671)
l3 = SphericalLens(85.87e-3, -646.31e-3, 7.05e-3, 60e-3, λ -> 1.4892)
l4 = Lens(SphericalSurface(-646.31e-3, 60e-3), SphericalSurface(23.51e-3, 40e-3), 1.9e-3, λ -> 1.7394)
translate3d!(l3, [0, thickness(l2), 0])
translate3d!(l4, [0, thickness(l2) + thickness(l3), 0])
l234 = TripletLens(l2, l3, l4)
l567 = SphericalTripletLens(Inf, 51.09e-3, -22.12e-3, -103.13e-3, 2.48e-3, 19.81e-3, 4.57e-3, 42e-3,
                            λ -> 1.5232, λ -> 1.6578, λ -> 1.5894)

# Calculate translation distances
l_234 = thickness(l1) + 0.38e-3
l_567 = l_234 + thickness(l234) + 13.0e-3 + 2.24e-3

# move elements into position
translate3d!(l234, [0, l_234, 0])
translate3d!(l567, [0, l_567, 0])

# spot detector in the paraxial image plane, back focal length is 44.902 mm
pd = Detector(50mm)
translate3d!(pd, [0, l_567 + thickness(l567) + 44.902mm, 0])

sonnar = ObjectGroup([l1, l234, l567])

system = StaticSystem([sonnar, pd])

##
λ = 587.6e-9 # m
fields = [
    (h = 0.0,    θ = 0.0,     z0 = 0.0,       color = :blue),
    (h = 21.6mm, θ = 12.0619, z0 = -25.319mm, color = :green),
    (h = 43.2mm, θ = 23.3057, z0 = -54.766mm, color = :red),
]

dir_theta(θ) = [0, cosd(θ), sind(θ)]
b1 = CollimatedSource([0, -50mm, fields[1].z0], dir_theta(fields[1].θ), 50mm, λ; num_rays=500)
b2 = CollimatedSource([0, -50mm, fields[2].z0], dir_theta(fields[2].θ), 40mm, λ; num_rays=500)

solve_system!(system, b1)
b1_spots = spot_diagram(pd)
empty!(pd)

solve_system!(system, b2)
b2_spots = spot_diagram(pd)
empty!(pd)

##
cview = [
  0.35614    0.934433   -2.08167e-17  -0.0304396
 -0.167338   0.0637774   0.983835      0.00380212
  0.919327  -0.350382    0.17908      -0.10951
  0.0        0.0         0.0           1.0
]

fig = Figure(size=(600, 600))
ax = LScene(fig[1,1:2], show_axis=false)

render!(ax, sonnar; transparency=true, alpha=0.1)
# render!(ax, stop)
# render!(ax, s7)
render!(ax, pd)

render!(ax, b1; alpha=0.15, render_every=5, color=:red)
render!(ax, b2; alpha=0.15, render_every=5, color=:green)

spot_ax1 = Axis(
    fig[2,1],
    aspect=1,
    xlabel="x [µm]",
    ylabel="z [µm]",
    yaxisposition = :left,
    title = "h = $(round(fields[1].h / mm, digits=1)) mm, θ = $(fields[1].θ)°",
)
scatter!(spot_ax1, b1_spots*1e6, markersize=4, color=:red)

spot_ax2 = Axis(
    fig[2,2],
    aspect=1,
    xlabel="x [µm]",
    ylabel="z [µm]",
    yaxisposition = :right,
    title = "h = $(round(fields[2].h / mm, digits=1)) mm, θ = $(round(fields[2].θ, digits=2))°",
)
scatter!(spot_ax2, b2_spots*1e6, markersize=4, color=:green)

display(fig)
set_view(ax, cview)

save("sonnar_spot_diagram.png", fig, px_per_unit=4, update=false)