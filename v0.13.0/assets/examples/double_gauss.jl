using GLMakie, CairoMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics
const mm = 1e-3

include(joinpath(@__DIR__, "..", "render_utils.jl"))

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

# generate render
cview = [
    0.484123   0.875     -1.11022e-16  -0.0272104
    -0.260453   0.144104   0.954672     -0.00137587
     0.835338  -0.462179   0.29766      -0.0651289
     0.0        0.0        0.0           1.0
]

fig = Figure(size=(600,380))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system)

λ = 486.0 # nm
zs = LinRange(-0.02, 0.02, 10)
for (i, z) in enumerate(zs)
    beam = Beam(Ray([0, -0.05, z], [0, 1, 0], λ))
    solve_system!(system, beam)
    render!(ax, beam, flen=0.06, show_pos=true)
end

set_view(ax, cview)
save("double_gauss.png", fig; px_per_unit=8, update = false)

## thin lens comparison
CairoMakie.activate!()

thin_lens = SphericalLens(100e-3, 100e-3, 0, 52.3e-3, λ -> 1.5)

tl_system = StaticSystem([thin_lens])

fig = Figure(size=(600,380))
aspect = (1,2,1)
limits = (-0.05, 0.05, -0.05, 0.15, -0.05, 0.05)
ax = Axis3(fig[1,1]; aspect, limits, azimuth=0, elevation=1e-3)

# hide decorations for vis. purposes
hidexdecorations!(ax)
hidezdecorations!(ax)

render!(ax, tl_system)

λ = 486.0 # nm
zs = LinRange(-0.02, 0.02, 10)
for (i, z) in enumerate(zs)
    beam = Beam(Ray([0, -0.05, z], [0, 1, 0], λ))
    solve_system!(tl_system, beam)
    render!(ax, beam, flen=0.2)
end

save("thin_lens_f100.png", fig, px_per_unit=4)