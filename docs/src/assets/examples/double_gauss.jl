using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

##
# define spherical lenses
l1 = SphericalLens(48.88e-3, 182.96e-3, 8.89e-3, 52.3e-3, λ -> 1.62286)
l23 = SphericalDoubletLens(36.92e-3, Inf, 23.06e-3, 15.11e-3, 2.31e-3, 45.11e-3, λ -> 1.58565, λ -> 1.67764)
l45 = SphericalDoubletLens(-23.91e-3, Inf, -36.92e-3, 1.92e-3, 7.77e-3, 40.01e-3, λ -> 1.57046, λ -> 1.64128)
l6 = SphericalLens(1063.24e-3, -48.88e-3, 6.73e-3, 45.11e-3, λ -> 1.62286)

# Calculate translation distances
l_23 = thickness(l1) + 0.38e-3
l_45 = l_23 + thickness(l23) + 9.14e-3 + 13.36e-3
l_6 = l_45 + thickness(l45) + 0.38e-3

# move elements into position
translate3d!(l23, [0, l_23, 0])
translate3d!(l45, [0, l_45, 0])
translate3d!(l6, [0, l_6, 0])

system = StaticSystem([l1, l23, l45, l6])

# generate render, same view for both lenses
cview = [
    0.484123   0.875     -1.11022e-16  -0.0322104
    -0.260453   0.144104   0.954672     -0.00137587
     0.835338  -0.462179   0.29766      -0.0951289
     0.0        0.0        0.0           1.0
]

fig = Figure(size=(600,380))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)

λ = 486e-9 # m
zs = LinRange(-0.02, 0.02, 10)
for (i, z) in enumerate(zs)
    beam = Beam(Ray([0, -0.05, z], [0, 1, 0], λ))
    solve_system!(system, beam)
    render!(ax, beam, flen=0.06, show_pos=true)
end

set_view(ax, cview)
save("double_gauss.png", fig; px_per_unit=8, update = false)

## Sonnar comparison
s1 = SphericalLens(69.21e-3, 433.84e-3, 9.33e-3, 70e-3, λ -> 1.671)
# front triplet: last surface only has a clear aperture of 40 mm -> assembled from individual lenses
s2 = SphericalLens(35.86e-3, 85.87e-3, 11.81e-3, 60e-3, λ -> 1.671)
s3 = SphericalLens(85.87e-3, -646.31e-3, 7.05e-3, 60e-3, λ -> 1.4892)
s4 = Lens(SphericalSurface(-646.31e-3, 60e-3), SphericalSurface(23.51e-3, 40e-3), 1.9e-3, λ -> 1.7394)
translate3d!(s3, [0, thickness(s2), 0])
translate3d!(s4, [0, thickness(s2) + thickness(s3), 0])
s234 = TripletLens(s2, s3, s4)
s567 = SphericalTripletLens(Inf, 51.09e-3, -22.12e-3, -103.13e-3, 2.48e-3, 19.81e-3, 4.57e-3, 42e-3,
                            λ -> 1.5232, λ -> 1.6578, λ -> 1.5894)

# Calculate translation distances
s_234 = thickness(s1) + 0.38e-3
s_567 = s_234 + thickness(s234) + 13.0e-3 + 2.24e-3

# move elements into position
translate3d!(s234, [0, s_234, 0])
translate3d!(s567, [0, s_567, 0])

sonnar = StaticSystem([s1, s234, s567])

fig = Figure(size=(600,380))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, sonnar)

# same relative filling of the aperture as above: F/1.5 instead of F/2
zs = LinRange(-0.0267, 0.0267, 10)
for (i, z) in enumerate(zs)
    beam = Beam(Ray([0, -0.05, z], [0, 1, 0], λ))
    solve_system!(sonnar, beam)
    render!(ax, beam, flen=0.045, show_pos=true)
end

set_view(ax, cview)
save("sonnar.png", fig; px_per_unit=8, update = false)
