using GLMakie, BeamletOptics

GLMakie.activate!()

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

##
λ = 1000nm
f = 5mm
D = 12mm

# on-axis paraboloid: vertex at the origin, focus at (0, -f, 0)
mirror = ParabolicMirror(f, 12.5mm, thickness=1mm)
# translate_to3d!(mirror, [0, 50mm, 0])

# spawn the beam between detector (y = -f) and mirror rim (y ≈ -1.95 mm), so the detector does not clip it
src = UniformDiscSource([0, -f / 1.1, 0], [0, 1, 0], D, λ; num_rays = 5000)
beam = CollimatedSource(
    [Beam(position(first(rays(b))), direction(first(rays(b))), λ, [1.0, 0, 0]) for b in BMO.beams(src)], D, [0, -f / 1.1, 0], [0, 1, 0])

pd = Detector(1mm)
translate_to3d!(pd, [0, -f, 0])

system = System([mirror, pd]);

solve_system!(system, beam)

##
NA = sin(2 * atan(D / (4f)))
R = 1.5λ / NA
xs, zs, E = electric_field(pd; n = 201, x_min = -R, x_max = R, z_min = -R, z_max = R)

I = intensity.(E)
I0 = maximum(I)
component(k) = map(e -> abs2(e[k]), E) ./ (2 * BMO.Z_vacuum * I0)

panels = [
    ("I (total)", I ./ I0),
    ("|E_x|² (along polarization)", component(1)),
    ("|E_y|² (longitudinal)", component(2)),
    ("|E_z|²", component(3)),
]

vector_psf_fig = Figure(size = (900, 1100))
for (k, (title, data)) in enumerate(panels)
    row, col = fldmod1(k, 2)
    ax = Axis(vector_psf_fig[row, col][1, 1]; title, xlabel = "x [µm]", ylabel = "z [µm]", aspect = DataAspect())
    hm = heatmap!(ax, xs * 1e6, zs * 1e6, data; colormap = :inferno, colorrange = (0, maximum(data)))
    Colorbar(vector_psf_fig[row, col][1, 2], hm)
end

c = (length(xs) + 1) ÷ 2
cut_ax = Axis(vector_psf_fig[3, 1:2]; title = "Cuts through the focus (parabolic mirror, NA = $(round(NA, digits = 2)))",
    xlabel = "position [µm]", ylabel = "normalized intensity")
lines!(cut_ax, xs * 1e6, I[:, c] ./ I0; label = "along x (polarization)")
lines!(cut_ax, zs * 1e6, I[c, :] ./ I0; label = "along z", linestyle = :dash)
axislegend(cut_ax)

save("psf_vector_showcase.png", vector_psf_fig, px_per_unit = 2)

##
fig = Figure(size = (600, 400))
ax = Axis3(fig[1,1]; azimuth=-0.2, elevation=1e-3)

hidedecorations!(ax)
hidespines!(ax)

render!(ax, mirror, transparency=true, alpha=0.5)
render!(ax, pd, color=:red)
render!(ax, beam; render_every=50, alpha=0.1, show_pos=false)

save("psf_vector_showcase2.png", fig, px_per_unit = 4)