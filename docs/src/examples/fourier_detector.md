# Fourier beamlet detectors

The [`FourierBeamletPropagator`](@ref) provides high-performance, non-uniform Fourier transform (NUFFT) based field synthesis for optical detectors in **BeamletOptics.jl**. It combines the efficiency of geometric ray and Gaussian beamlet tracing with wave-optical diffraction and interference effects.

## Physical principles and architecture

When rays or beamlets intersect an [`AbstractFourierDetector`](@ref) surface:

1. **Discrete Fourier k-points:** Each incident ray or beamlet is converted into a [`FourierKPoint`](@ref) containing its 3D wavevector \$\mathbf{k}\$, optical frequency \$\omega\$, complex electric field amplitude \$\mathbf{E}_0\$, and optical path length (OPL) phase anchor.
2. **Frequency-space field synthesis:** Using `synthesize_field` or [`calculate_observables`](@ref), the coherent electromagnetic field is evaluated over an arbitrary 1D, 2D, or 3D [`SpatialGrid`](@ref) via accelerated non-uniform FFTs (`NonuniformFFTs.jl`) or direct phasor accumulation:
   ```math
   \mathbf{E}(\mathbf{r}) = \sum_j \mathbf{E}_j \exp\bigl(i \mathbf{k}_j \cdot \mathbf{r} + i \phi_j\bigr)
   ```
3. **Polarimetric observables:** The field synthesis computes full 3D vectorial observables in a [`DetectorObservables`](@ref) struct:
   - `intensity`: Poynting flux \$I(\mathbf{r}) \ge 0\$ in \$\text{W/m}^2\$.
   - `phase`: Wavefront phase \$\phi(\mathbf{r}) \in [-\pi, \pi]\$ in radians.
   - `stokes`: Full 4-vector Stokes polarization parameters \$(S_0, S_1, S_2, S_3)\$.
   - `power`: Total integrated optical power in Watts (\$P = \iint I \, dA\$).

The examples below illustrate three standard configurations: diffraction-limited focal spot estimation (Airy disc), two-beam interferometric fringe formation, and 1D anamorphic focusing with a cylindrical lens.

---

## 1. Spherical singlet focus and Airy disc

When a collimated, circular beam passes through a spherical focusing lens, diffractive aperture truncation at the pupil forms an Airy diffraction pattern at the focal plane. Under paraxial approximations, the intensity profile follows:

```math
I(r) = I_0 \left[ \frac{2 J_1(\pi D r / (\lambda f))}{\pi D r / (\lambda f)} \right]^2
```

with first dark ring radius:

```math
r_{\text{Airy}} \approx 1.22 \frac{\lambda f}{D}
```

The [`FourierBeamletPropagator`](@ref) collects the refracted rays with their exact accumulated optical path lengths and synthesizes the focal diffraction spot.

```julia
using BeamletOptics
using LinearAlgebra

# 1. System parameters
R = 103.36e-3       # Front radius of curvature (m)
n = 1.5168          # Refractive index (N-BK7)
d_lens = 25.4e-3    # Lens aperture diameter (m)
l_lens = 1.0e-3     # Center thickness (m)
D_beam = 10.0e-3    # Beam diameter (m)
λ = 632.8e-9        # HeNe laser wavelength (m)
f = 200.0e-3        # Effective focal length (m)
y_foc = 200.13e-3   # Focal plane position along optical axis +y (m)

# 2. Components
src = UniformDiscSource([0.0, -10.0e-3, 0.0], [0.0, 1.0, 0.0], D_beam, λ; num_rays = 2000)
lens = SphericalLens(R, Inf, l_lens, d_lens, x -> n)

# Planar Fourier detector positioned at the focal plane
det = FourierBeamletPropagator(100.0e-6; stop = true, is_planar = true)
translate3d!(det, [0.0, y_foc, 0.0])

# 3. Ray tracing and solving
sys = System([lens, det])
solve_system!(sys, src)

# 4. Observables evaluation on a high-resolution grid
xs = range(-50.0e-6, 50.0e-6, length = 129)
zs = range(-50.0e-6, 50.0e-6, length = 129)
obs = calculate_observables(det, (xs, zs))

# Extract intensity and verify central peak and Airy minimum
I = obs.intensity
I_max = maximum(I)
mid = (129 + 1) ÷ 2

radial_profile = I[mid:end, mid]
min_idx = argmin(radial_profile[1:30])
r_min = (min_idx - 1) * step(xs)

println("Peak intensity: ", I[mid, mid])
println("First dark ring radius: ", r_min * 1e6, " µm (theoretical: 15.44 µm)")
```

---

## 2. Mach-Zehnder tilt fringes

In a Mach-Zehnder interferometer, a coherent beam is split into two separate paths by a 50:50 beamsplitter and recombined at a second beamsplitter. When one folding mirror is tilted by a small angle \$\Delta\theta\$, the two recombining wavefronts acquire a relative inclination, producing straight, equally spaced interference fringes:

```math
\Lambda = \frac{\lambda}{\Delta\theta}
```

The [`FourierBeamletPropagator`](@ref) accumulates the coherent fields from both optical arms without mutual cross-talk and reconstructs the interference fringes with full visibility.

```julia
using BeamletOptics
using LinearAlgebra

# 1. Geometry setup (2-inch baseline)
inch = BeamletOptics.inch
m1 = SquarePlanoMirror2D(inch)
m2 = SquarePlanoMirror2D(inch)
b1 = ThinBeamsplitter(inch, reflectance = 0.5)
b2 = ThinBeamsplitter(inch, reflectance = 0.5)

translate3d!(b1, [0 * inch, 0 * inch, 0.0])
translate3d!(b2, [2 * inch, 2 * inch, 0.0])
translate3d!(m1, [0 * inch, 2 * inch, 0.0])
translate3d!(m2, [2 * inch, 0 * inch, 0.0])

# 45-degree alignments
zrotate3d!(b1, deg2rad(360 - 135))
zrotate3d!(b2, deg2rad(45))
zrotate3d!(m1, deg2rad(360 - 135))

# Introduce a small tilt in mirror m2 around z-axis
tilt_angle = deg2rad(0.02)
zrotate3d!(m2, deg2rad(45) + tilt_angle)

# 2. Place detector at output port along +y
det = FourierBeamletPropagator(10.0e-3; stop = true, is_planar = true)
translate3d!(det, [2 * inch, 3 * inch, 0.0])

# 3. Trace polarized beam through interferometer
λ = 632.8e-9
ray = PolarizedRay([0.0, -0.05, 0.0], [0.0, 1.0, 0.0], λ, [0.0, 0.0, 1.0])
beam = Beam(ray)

sys = System([b1, m1, m2, b2, det])
solve_system!(sys, beam)

# 4. Evaluate fringe pattern across detector aperture
xs = range(-2.0e-3, 2.0e-3, length = 101)
zs = range(-2.0e-3, 2.0e-3, length = 101)
obs = calculate_observables(det, (xs, zs))

# Fringe visibility
I_max = maximum(obs.intensity)
I_min = minimum(obs.intensity)
visibility = (I_max - I_min) / (I_max + I_min)

println("Interferometer fringe visibility: ", round(visibility, digits = 4))
```

---

## 3. Cylindrical lens astigmatic focus

A cylindrical lens features curvature in one dimension and zero optical power in the orthogonal dimension. An incident collimated circular beam is focused into an anamorphic focal line rather than a point, introducing strong 1D astigmatism.

By evaluating the synthesized field on a [`SpatialGrid`](@ref), the focal line characteristics—tight focus along the curved coordinate (\$z\$) and broad, unconstrained extent along the flat coordinate (\$x\$)—can be quantitatively extracted.

```julia
using BeamletOptics
using LinearAlgebra

# 1. Cylindrical lens parameters
r_cyl = 25.0e-3     # Cylinder radius of curvature (m)
d_cyl = 25.0e-3     # Width (m)
h_cyl = 25.0e-3     # Height (m)
ct_cyl = 5.0e-3     # Center thickness (m)
n_cyl = 1.5168      # N-BK7 refractive index
λ = 632.8e-9

# Focal length in the curved dimension: f ≈ r / (n - 1)
f_cyl = r_cyl / (n_cyl - 1.0) # ≈ 48.4 mm

cyl_surface = CylindricalSurface(r_cyl, d_cyl, h_cyl)
cyl_lens = Lens(cyl_surface, ct_cyl, x -> n_cyl)

# 2. Position planar detector at the focal distance
det = FourierBeamletPropagator(20.0e-3; stop = true, is_planar = true)
translate3d!(det, [0.0, f_cyl + ct_cyl, 0.0])

# 3. Collimated disc source
src = UniformDiscSource([0.0, -10.0e-3, 0.0], [0.0, 1.0, 0.0], 10.0e-3, λ; num_rays = 1000)

sys = System([cyl_lens, det])
solve_system!(sys, src)

# 4. Compute field observables across (x, z) plane
xs = range(-5.0e-3, 5.0e-3, length = 129)
zs = range(-2.0e-3, 2.0e-3, length = 129)
obs = calculate_observables(det, (xs, zs))

# Compare line widths: narrow waist along focused z, wide extent along unfocused x
println("Line focus maximum intensity: ", maximum(obs.intensity))
println("Total collected power: ", obs.power, " W")
```

---

## Resetting detectors for multiple runs

As with standard [`Detector`](@ref) instances, an [`AbstractFourierDetector`](@ref) accumulates hits across multiple solve steps or beams until explicitly reset:

```julia
# Reset captured hits between parameter sweeps
empty!(det)
```
