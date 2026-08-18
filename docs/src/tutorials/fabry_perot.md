# Fabry-Pérot Cavity Resonator

This tutorial demonstrates how to simulate a Fabry-Pérot cavity resonator and model coherent multiple-reflection interference using the ray-splitting features of `BeamletOptics.jl`.

We will explore:
- Setting up a parallel-mirror Fabry-Pérot etalon.
- Tracing an `AstigmaticGaussianBeamlet` with active ray-splitting.
- Pruning infinite round trips using the `power_cutoff` option.
- Simulating the wavelength transmission spectrum (Airy distribution).
- Introducing a small wedge to observe spatial interference fringes (Fizeau-like fringes).

---

## Cavity Resonator Concept

A Fabry-Pérot resonator consists of a transparent plate (etalon) or cavity bounded by two parallel, partially reflective surfaces. 

When a light beam is incident on the cavity:
- A portion of the light enters the cavity through the front face.
- Inside the cavity, the light undergoes multiple internal reflections between the front and back surfaces.
- At each bounce, a fraction of the light transmits through the back face.
- These multiple transmitted beams exit parallel to each other and interfere coherently.

If the round-trip path length is an integer multiple of the laser wavelength, the transmitted beams interfere constructively, resulting in maximum transmission (resonance). Otherwise, they interfere destructively, causing high reflectance.

---

## Parallel Mirror Resonance Scan

First, we set up a parallel-mirror cavity (etalon) using a flat glass substrate (modelled as a plano-planar `SphericalLens` with infinite radii of curvature) coated on both sides with a beam-splitting coating.

```@example fabry_perot
using BeamletOptics
using LinearAlgebra
using CairoMakie


const mm = 1e-3
const μm = 1e-6

# Define cavity parameters
d = 10.0μm          # Center spacing
n_glass = 1.5       # Refractive index of substrate
R_power = 0.8       # Power reflectance (80%)
T_power = 1.0 - R_power

r = sqrt(R_power)
t = sqrt(T_power)

# Define a beam-splitting coating
splitting_coating = SimpleBeamsplitterCoating(r, r, t, t)

# Construct flat glass etalon (SphericalLens of thickness d)
etalon_base = SphericalLens(Inf, Inf, d, 2.0mm, λ -> n_glass)
etalon = etalon_base |> with_coatings(front = splitting_coating, back = splitting_coating)

# Place a photodetector behind the etalon
pd = Detector(1.0mm)
translate3d!(pd, [0.0, 1.0mm, 0.0])

# Assemble the system
system = System([etalon, pd])
nothing # hide
```

Since the reflections continue infinitely, a simple ray trace would lead to infinite branching. To keep the simulation finite, we use the `power_cutoff = 1e-4` argument in `solve_system!`. Once the power weight of a split ray path falls below $0.01\%$, tracing on that branch is stopped and the sub-branch is dropped.

We will scan the wavelength of the incident astigmatic Gaussian beamlet from 620 nm to 630 nm to observe the Airy resonance peak.

```@example fabry_perot
lambdas = LinRange(620e-9, 630e-9, 101)
powers = Float64[]

for λ in lambdas
    empty!(pd)
    # Start AstigmaticGaussianBeamlet propagating along +y with explicit initial power P0 = 1.0 W
    beam = AstigmaticGaussianBeamlet([0.0, -1.0mm, 0.0], [0.0, 1.0, 0.0], λ, 100.0μm, P0 = 1.0)
    solve_system!(system, beam; power_cutoff = 1e-4, retrace = false)
    push!(powers, optical_power(pd; n=20))
end

# Find resonance peaks
fig_airy = Figure(size = (600, 350))
ax_airy = Axis(fig_airy[1, 1], xlabel = "Wavelength [nm]", ylabel = "Transmitted Power [W]", title = "Airy Resonance Spectrum")
lines!(ax_airy, lambdas ./ 1e-9, powers, color = :red, linewidth = 2)
vlines!(ax_airy, 625.0, color = :blue, linestyle = :dash) # Theoretical resonance peak (m = 48)

save("fabry_perot_airy.png", fig_airy) # hide
nothing # hide
```

![Airy Resonance Peak](fabry_perot_airy.png)

At $\lambda = 625$ nm, the optical path length of a round trip inside the glass is $2 \cdot n \cdot d = 2 \cdot 1.5 \cdot 10\,\mu\text{m} = 30\,\mu\text{m}$, which is exactly $48$ times the wavelength. Thus, the multiple internal reflections interfere constructively, producing a sharp transmission peak of over $90\%$ of the initial laser power.

---

## Wedged Cavity and Spatial Fizeau Fringes

Next, we introduce a slight wedge to the cavity. Instead of parallel mirrors, we tilt the mirrors relative to each other by $0.1^\circ$. This creates a spatially varying cavity thickness along the x-direction:
$$d(x) = d_0 + x \cdot \tan\alpha$$

This thickness variation translates to a spatially varying phase shift across the profile of the beam, producing spatial interference fringes (Fizeau fringes) on the detector face.

### Physical Astigmatism from Mirror Tilt

Introducing a wedge angle means the beam's reflections off the cavity mirrors are no longer at normal incidence. Non-normal reflection off curved or flat boundaries introduces astigmatism to the wavefront (the tangential and sagittal curvatures propagate differently). 

To model this physical phenomenon accurately, we must use an `AstigmaticGaussianBeamlet` (AGB) rather than a circular `GaussianBeamlet`. The AGB dynamically updates its full 3D complex curvature matrix at each tilted surface crossing, capturing the realistic spatial profiles of the interfering beams.

To build the wedged cavity, we construct the setup using two separate flat `SphericalLens` objects with coatings so we can orient them independently:

```@example fabry_perot
# Wedge angle
α = deg2rad(0.1)

# Mirror 1 (Entrance) - Flat CoatedLens tilted by +α/2
m1_base = SphericalLens(Inf, Inf, 0.5mm, 5.0mm, λ -> n_glass)
mirror1 = m1_base |> with_coatings(back = splitting_coating, front = SimpleARCoating(0.0))
translate3d!(mirror1, [0.0, -0.25mm, 0.0]) # Front face near y = -0.5mm, back face at y = 0
zrotate3d!(mirror1, α/2)

# Mirror 2 (Exit) - Flat CoatedLens tilted by -α/2
m2_base = SphericalLens(Inf, Inf, 0.5mm, 5.0mm, λ -> n_glass)
mirror2 = m2_base |> with_coatings(front = splitting_coating, back = SimpleARCoating(0.0))
translate3d!(mirror2, [0.0, 20.0μm + 0.25mm, 0.0]) # Front face at y = 20μm
zrotate3d!(mirror2, -α/2)

# Photodetector placed behind the cavity
pd_wedge = Detector(1.5mm)
translate3d!(pd_wedge, [0.0, 2.0mm, 0.0])

# Assemble the wedged system
system_wedge = System([mirror1, mirror2, pd_wedge])
```

Now, we launch a wide astigmatic Gaussian beamlet ($\lambda = 632.8$ nm, waist $w_0 = 300\,\mu$m, explicitly setting initial power $P_0 = 1.0$ W) and capture the high-resolution 2D spatial intensity map on the detector.

```@example fabry_perot
# Start beam and solve the system
beam_wedge = AstigmaticGaussianBeamlet([0.0, -1.0mm, 0.0], [0.0, 1.0, 0.0], 632.8e-9, 300.0μm, P0 = 1.0)
solve_system!(system_wedge, beam_wedge; power_cutoff = 1e-4, retrace = false)

# Compute high-resolution 2D intensity grid
x_grid, z_grid, I_map = intensity(pd_wedge, n=500)

# Plot Fizeau fringes
fig_fizeau = Figure(size = (600, 500))
ax_fizeau = Axis(fig_fizeau[1, 1], xlabel = "x [mm]", ylabel = "z [mm]", title = "Fizeau Spatial Interference Fringes", aspect = 1)
heatmap!(ax_fizeau, x_grid ./ mm, z_grid ./ mm, I_map, colormap = :viridis)

save("fabry_perot_fizeau.png", fig_fizeau) # hide
nothing # hide
```

![Fizeau Fringes](fabry_perot_fizeau.png)

The heatmap displays parallel vertical fringes, demonstrating how the phase differences vary linearly across the beam profile due to the wedge geometry. The fringe spacing matches the analytical expectation $\Delta x = \lambda / (2 \cdot \tan\alpha) \approx 0.18$ mm perfectly.
