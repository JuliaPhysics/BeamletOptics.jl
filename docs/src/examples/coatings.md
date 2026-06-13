# Optical Coatings Examples

This page provides practical examples demonstrating how to use coatings in `BeamletOptics.jl`.

---

## Quarter-Wave Anti-Reflective Coating

A standard method to minimize reflection at a dielectric boundary is a single quarter-wave thin-film layer. For normal incidence, the optimal refractive index of the coating layer $n_c$ between medium $n_1$ and $n_2$ is:

```math
n_c = \sqrt{n_1 n_2}
```

The optimal physical thickness $d$ of the layer is:

```math
d = \frac{\lambda_0}{4 n_c}
```

where $\lambda_0$ is the design wavelength.

The following example compares the transmission of a single ray through an uncoated lens versus an AR-coated lens:

```@example coatings_showcase
using BeamletOptics
const BMO = BeamletOptics

# Setup parameters
λ = 1000e-9        # 1000 nm design wavelength
n_air = 1.0
n_glass = 1.5      # substrate index

# Calculate optimal AR coating properties
n_c = sqrt(n_air * n_glass)   # ≈ 1.2247
d_c = λ / (4 * n_c)           # ≈ 204 nm

# Define the coating model
ar_coating = ThinFilmCoating(n_c, d_c)

# Create uncoated and coated plano-planar lenses
glass_index = λ -> n_glass
lens_uncoated = SphericalLens(Inf, Inf, 5e-3, 10e-3, glass_index)
lens_coated = lens_uncoated |> with_coatings(front = ar_coating)

# Trace a polarized ray through both systems
ray_in_uncoated = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])
ray_in_coated = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], λ, [1.0, 0.0, 0.0])

sys_uncoated = System([lens_uncoated])
beam_uncoated = Beam(ray_in_uncoated)
solve_system!(sys_uncoated, beam_uncoated; retrace = false)

sys_coated = System([lens_coated])
beam_coated = Beam(ray_in_coated)
solve_system!(sys_coated, beam_coated; retrace = false)

# Compare transmission at the front surface (Ray 2 is inside the lens)
ray_inside_uncoated = rays(beam_uncoated)[2]
ray_inside_coated = rays(beam_coated)[2]

# Uncoated amplitude transmission: 2 * n1 / (n1 + n2) = 2 * 1 / (1 + 1.5) = 0.8
println("Uncoated inside field amplitude: ", abs(BMO.polarization(ray_inside_uncoated)[1]))

# Coated amplitude transmission: 100% power transmission -> amplitude factor is sqrt(n1 / n2) ≈ 0.8165
println("AR-coated inside field amplitude: ", abs(BMO.polarization(ray_inside_coated)[1]))
```

---

## Polarizing Cube Beamsplitter

Prism cube beamsplitters use a split-behavior coating on the hypotenuse interface of two cemented right-angle prisms. Here, we model a 50:50 splitter:

```@example coatings_showcase
# Define leg length and substrate index
leg_length = 10e-3
glass_idx = DiscreteRefractiveIndex([λ], [1.5])

# Construct a CubeBeamsplitter with 50% reflectance
cube_splitter = CubeBeamsplitter(leg_length, glass_idx; reflectance = 0.5)

# Trace an Astigmatic Gaussian Beamlet at normal incidence to the front face
agb = AstigmaticGaussianBeamlet(
    [0.0, -15e-3, 0.0], [0.0, 1.0, 0.0], λ, 1e-3;
    E0 = [1.0, 0.0, 0.0], support = [0.0, 0.0, 1.0]
)

system_bs = System([cube_splitter])
solve_system!(system_bs, agb; retrace = false)

# The beamsplitter splits the beamlet into a transmitted and reflected child
println("Number of split child beams: ", length(agb.children))

t_beam = agb.children[1]
r_beam = agb.children[2]

# Show resulting intensities/amplitudes
using LinearAlgebra
println("Transmitted beam field amplitude: ", norm(BMO.electric_field(t_beam)))
println("Reflected beam field amplitude: ", norm(BMO.electric_field(r_beam)))
```

---

## Multi-layer Dielectric bandpass Filter

Multi-layer thin-film coatings are used to make highly wavelength-selective optical filters (such as bandpass filters). This is done by stacking alternating layers of high-index (e.g. $TiO_2$, $n_H = 2.4$) and low-index (e.g. $SiO_2$, $n_L = 1.45$) materials.

The following example designs a multi-layer stack of 9 alternating quarter-wave layers centered at $\lambda_0 = 800~\text{nm}$ and plots its reflectance spectrum:

```@example coatings_showcase
using CairoMakie

# Design parameters
λ0 = 800e-9
nH = 2.4
nL = 1.45

# Stacking alternating H and L layers
# HLHLHLHLH
ns = [nH, nL, nH, nL, nH, nL, nH, nL, nH]
ds = [λ0 / (4 * n) for n in ns]

# Construct the thin-film coating model
filter_coating = ThinFilmCoating(ns, ds; behavior = BMO.Reflective())

# Evaluate reflectance over a spectrum range
wavelengths = range(400e-9, 1200e-9, length=200)
reflectances = Float64[]

for wl in wavelengths
    rs, rp, ts, tp = BMO.fresnel_coefficients(filter_coating, 0.0, wl, 1.0, 1.5)
    push!(reflectances, abs(rs)^2)
end

# Plot the spectrum
fig = Figure(size=(600, 350))
ax = Axis(fig[1, 1], xlabel="Wavelength (nm)", ylabel="Reflectance", title="9-Layer Dielectric Filter Spectrum")
lines!(ax, wavelengths .* 1e9, reflectances, linewidth=2, color=:dodgerblue)
save("dielectric_filter_spectrum.png", fig); nothing # hide
```

![Dielectric filter spectrum](dielectric_filter_spectrum.png)

---

## Frustrated Total Internal Reflection (FTIR) / Optical Tunneling

When a ray strikes a glass-air boundary at an angle of incidence $\theta_i$ greater than the critical angle $\theta_c = \arcsin(1/n_{\text{glass}})$, total internal reflection (TIR) occurs. However, if a second glass block is placed extremely close to the first (within a fraction of a wavelength), the evanescent wave in the air gap is "frustrated". Light waves tunnel through the thin barrier, transmitting power into the second glass block.

We can model the thin air gap as a `ThinFilmCoating` with a splitting behavior. Below, we calculate and plot the transmission coefficient $T_s = |t_s|^2$ as a function of the air gap thickness $d$:

```@example coatings_showcase
# Setup physical parameters
λ = 632.8e-9        # 632.8 nm (He-Ne laser)
n_glass = 1.5
n_gap = 1.0         # Air gap
θi = deg2rad(45.0)  # 45° angle (> critical angle of ≈ 41.8°)

# Evaluate transmission over a range of gap thicknesses
gap_thicknesses = range(0.0, 800e-9, length=200)
transmissions = Float64[]

for d in gap_thicknesses
    # Model the air gap of thickness d
    gap_coating = ThinFilmCoating([n_gap], [d]; behavior = BMO.Splitting())
    rs, rp, ts, tp = BMO.fresnel_coefficients(gap_coating, θi, λ, n_glass, n_glass)
    
    # Power transmission transmissivity is |ts|^2 (since indices are matched)
    push!(transmissions, abs(ts)^2)
end

# Plot transmission vs air gap thickness
fig2 = Figure(size=(600, 350))
ax2 = Axis(fig2[1, 1], xlabel="Air Gap Thickness (nm)", ylabel="Power Transmission (T)", title="Frustrated Total Internal Reflection (Optical Tunneling)")
lines!(ax2, gap_thicknesses .* 1e9, transmissions, linewidth=2, color=:crimson)
save("ftir_transmission_curve.png", fig2); nothing # hide
```

![Frustrated Total Internal Reflection](ftir_transmission_curve.png)

