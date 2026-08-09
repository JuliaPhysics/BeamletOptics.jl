# Wedge interference filter light trap (WOLIT)

This example demonstrates how to implement a custom angle-dependent thin film interference filter coating and simulate a He-Ne laser wedge light trap (Wedge Optical Light Trap / WOLIT) as described by Ewers & Lorbeer (2019) [Ewers:2019](@cite).

---

## Wedge Light Trap Concept

A wedge light trap consists of two surfaces set at a small angle relative to each other (a wedge):
* The top surface is a highly-reflective (HR) mirror.
* The bottom surface is coated with a custom bandpass or edge filter.

The edge filter's cut-off wavelength shifts with the angle of incidence (AOI) according to the formula:

```math
\lambda_{\text{edge}}(\theta) = \lambda_0 \sqrt{1 - \frac{\sin^2(\theta_{\text{air}})}{n_{\text{eff}}^2}}
```

where $\lambda_0$ is the cut-off wavelength at normal incidence, $\theta_{\text{air}}$ is the angle of incidence in air, and $n_{\text{eff}}$ is the effective refractive index of the filter cavity.

A laser ray enters the wedge at a high angle of incidence where the filter is transmissive. As the ray bounces between the mirror and the filter, the wedge geometry reduces the angle of incidence with each bounce. Eventually, the angle of incidence falls below the critical threshold where the filter becomes highly reflective, trapping the ray. The ray travels towards the wedge apex, reverses direction, and travels back towards the entrance where the angle of incidence increases again until it exits the cavity.

---

## Custom Coating Definition

First, we define a custom coating type that wraps a transmissive model (AR coating) and a reflective model (HR coating), switching between them based on the angle-dependent cut-off wavelength.

```@example wolit_trap
using BeamletOptics
using LinearAlgebra
using CairoMakie


const mm = 1e-3

# Custom Coating definition representing the WOLIT filter
struct WolitFilterCoating
    lambda_0::Float64
    n_eff::Float64
    laser_lambda::Float64
    transmissive_model::SimpleARCoating
    reflective_model::SimpleHRCoating
end

# Constructor
function WolitFilterCoating(lambda_0, n_eff, laser_lambda)
    return WolitFilterCoating(
        lambda_0, n_eff, laser_lambda,
        SimpleARCoating(0.02),
        SimpleHRCoating(0.999)
    )
end

# Overload coating_behavior dynamically based on the ray's angle of incidence
function BeamletOptics.coating_behavior(coating::WolitFilterCoating, ray)
    int = BeamletOptics.intersection(ray)
    normal = BeamletOptics.normal3d(int)
    n_incident = BeamletOptics.refractive_index(ray)
    cos_θ = abs(dot(BeamletOptics.direction(ray), normal))
    sin_θ = sqrt(max(0.0, 1.0 - cos_θ^2))
    sin_θ_air = n_incident * sin_θ
    sin_θ_air = clamp(sin_θ_air, 0.0, 1.0)
    λ_edge = coating.lambda_0 * sqrt(1.0 - sin_θ_air^2 / coating.n_eff^2)
    
    if coating.laser_lambda > λ_edge
        return BeamletOptics.Transmissive()
    else
        return BeamletOptics.Reflective()
    end
end

# Overload get_jones_matrix
function BeamletOptics.get_jones_matrix(coating::WolitFilterCoating, θi, λ, n1, n2, is_reflected; from_front::Bool = true)
    # The filter formula uses the angle of incidence in air
    sin_θ_air = n1 * sin(θi)
    sin_θ_air = clamp(sin_θ_air, 0.0, 1.0)
    λ_edge = coating.lambda_0 * sqrt(1.0 - sin_θ_air^2 / coating.n_eff^2)
    
    if coating.laser_lambda > λ_edge
        # Transmissive state
        return BeamletOptics.get_jones_matrix(coating.transmissive_model, θi, λ, n1, n2, is_reflected; from_front = from_front)
    else
        # Reflective state
        return BeamletOptics.get_jones_matrix(coating.reflective_model, θi, λ, n1, n2, is_reflected; from_front = from_front)
    end
end
```

---

## Simulating the Wedge Cavity

We assemble the optical system containing the flat top mirror and the bottom wedge lens. We then trace a ray entering at an angle of 26.5 degrees.

```@example wolit_trap
# Setup wavelength and coating properties
λ = 632.8e-9         # He-Ne laser wavelength (red light)
λ_cut = 650e-9       # Filter cut-off at normal incidence (650 nm)
n_eff = 1.85         # Effective refractive index of the filter cavity
α = deg2rad(1.2)     # Wedge angle: 1.2 degrees
d = 1.0mm            # Wedge thickness at origin (x=0)

# Define custom WolitFilterCoating
wolit_coating = WolitFilterCoating(λ_cut, n_eff, λ)

# Mirror 1: highly reflective mirror at y = 1.0mm, tilted by -α/2
m1_base = RectangularPlanoMirror(15mm, 2mm, 0.25mm)
mirror1 = with_coatings(m1_base; front = SimpleHRCoating(0.999))
translate3d!(mirror1, [0.0, d, 0.0])
zrotate3d!(mirror1, -α/2)

# Mirror 2: custom coated flat glass lens substrate (filter) at y = 0, size 20mm
# The bottom face is AR coated and the top face is the WolitFilterCoating
#m2_base = Lens(RectangularFlatSurface(15mm), RectangularFlatSurface(15mm), 1mm, λ -> 1.5)
m2_base = Lens(CylindricalSurface(Inf, 2mm, 15mm), CylindricalSurface(Inf, 2mm, 15mm), 0.5mm, λ -> 1.5)
coated_lens2 = CoatedLens(m2_base; front = SimpleARCoating(0.0), back = wolit_coating)
translate3d!(coated_lens2, [0.0, -1.0mm, 0.0])
zrotate3d!(coated_lens2, α/2)

# Assemble the system
system = System([mirror1, coated_lens2])

# Define and trace the incoming ray
θ_in = deg2rad(26.5)
pos_in = [-5.5mm, -2.5mm, 0.0]
dir_in = [sin(θ_in), cos(θ_in), 0.0]
pol_in = [0.0, 0.0, 1.0] # s-polarized, orthogonal to x-y plane

ray_in = PolarizedRay(pos_in, dir_in, λ, pol_in)
beam = Beam(ray_in)

# Trace the ray through the system
solve_system!(system, beam; r_max = 200, retrace = false)
all_rays = rays(beam)
n_rays = length(all_rays)
println("Total segments traced: $n_rays")
```

---

## 3D Visualization

Using CairoMakie, we render the layout of the wedge trap and the path of the multi-bounce ray inside the cavity.

```@example wolit_trap
fig = Figure(size = (800, 600))
ax = LScene(fig[1, 1])
ax.show_axis[] = false

# Render the optical system and the traced beam
render!(ax, system)
render!(ax, beam, linewidth = 2.5, color = :red, flen=1mm)

c_view = [ 0.0296277   0.999555   -0.00357691  -0.000779786 # hide
 -0.558722    0.0195281   0.829125     0.00190157 # hide
  0.828825   -0.0225666   0.559052    -0.0131178 # hide
  0.0         0.0         0.0          1.0] # hide

ax.scene.camera.view[] = c_view # hide 
cam = ax.scene.camera_controls # hide
update_cam!(ax.scene, cam)     # hide


save("wolit_trap.png", fig) # hide
nothing # hide
```

![WOLIT Wedge Trap](wolit_trap.png)
