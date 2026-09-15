using CairoMakie, BeamletOptics, LinearAlgebra

CairoMakie.activate!()

const BMO = BeamletOptics

λ = 1e-6

# Hit of a ray converging on `focus`, with zero phase at the focus
function focus_hit(dir, E0; focus = zeros(3))
    ray = PolarizedRay(focus .- λ .* dir, dir, λ, E0)
    BMO.intersection!(ray, BMO.Intersection(λ, BMO.Point3(-dir)))
    return BMO.PolarizedRayHit(ray, BMO.optical_path_length(ray))
end

function aplanatic_focus(NA; N = 60, pol = [1.0, 0, 0])
    pd = Detector(1.0)
    for i in 1:N, j in 1:(4N)
        α = asin(NA) * (i - 0.5) / N
        φ = 2π * (j - 0.5) / (4N)
        dir = [sin(α) * cos(φ), cos(α), sin(α) * sin(φ)]
        e_r = [cos(φ), 0, sin(φ)]
        e_φ = [-sin(φ), 0, cos(φ)]
        e_ρ = [cos(α) * cos(φ), -sin(α), cos(α) * sin(φ)]
        E0 = sqrt(cos(α)) * sin(α) .* (dot(pol, e_r) .* e_ρ .+ dot(pol, e_φ) .* e_φ)
        push!(pd, focus_hit(dir, E0))
    end
    return pd
end

NA = 0.9
R = 1.5λ / NA
pd = aplanatic_focus(NA)
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
cut_ax = Axis(vector_psf_fig[3, 1:2]; title = "Cuts through the focus (NA = $NA)",
    xlabel = "position [µm]", ylabel = "normalized intensity")
lines!(cut_ax, xs * 1e6, I[:, c] ./ I0; label = "along x (polarization)")
lines!(cut_ax, zs * 1e6, I[c, :] ./ I0; label = "along z", linestyle = :dash)
axislegend(cut_ax)

save("psf_vector_showcase.png", vector_psf_fig, px_per_unit = 2)
