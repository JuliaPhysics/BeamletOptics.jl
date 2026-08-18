# Materials, Absorption & Gain

In `BeamletOptics.jl`, the optical properties of bulk media (such as glass lenses, crystals, gases, or absorbing filters) are modeled through the [`RefractiveIndex`](@ref) interface. Bulk media support both **real** (lossless) and **complex** (absorbing or amplifying) refractive indices.

---

## Complex Refractive Index & Physical Conventions

A bulk medium is characterized by its complex refractive index:

```math
\tilde{n}(\lambda) = n(\lambda) + i \kappa(\lambda)
```

where:
* ``n(\lambda) = \text{Re}(\tilde{n})`` is the **real refractive index**, governing phase velocity and Snellius refraction.
* ``\kappa(\lambda) = \text{Im}(\tilde{n})`` is the **extinction coefficient**, governing attenuation or amplification.

### Sign Conventions

`BeamletOptics.jl` follows standard optical conventions (where plane waves propagate as ``e^{i(kz - \omega t)}``):

* **Bulk Absorption (``\kappa > 0``):**
  The linear power absorption coefficient ``\alpha(\lambda)`` in $[1/\text{m}]$ (Lambert-Beer law) is given by:
  ```math
  \alpha(\lambda) = \frac{4\pi \kappa(\lambda)}{\lambda}
  ```
  * **Power Attenuation ([`Ray`](@ref)):**
    ```math
    w(L) = w_0 \cdot e^{-\alpha L} = w_0 \cdot e^{-\frac{4\pi \kappa}{\lambda} L}
    ```
  * **Field Amplitude Attenuation ([`PolarizedRay`](@ref), [`GaussianBeamlet`](@ref), [`AstigmaticGaussianBeamlet`](@ref)):**
    ```math
    E(L) = E_0 \cdot e^{-\frac{\alpha}{2} L} = E_0 \cdot e^{-\frac{2\pi \kappa}{\lambda} L}
    ```

* **Small-Signal Gain (``\kappa < 0``):**
  For active laser media with small-signal power gain coefficient ``g(\lambda)`` in $[1/\text{m}]$:
  ```math
  g(\lambda) = -\alpha(\lambda) = -\frac{4\pi \kappa(\lambda)}{\lambda}
  ```
  * **Power Amplification:** ``w(L) = w_0 \cdot e^{+g L}``
  * **Field Amplification:** ``E(L) = E_0 \cdot e^{+\frac{g}{2} L}``

* **Snell's Law & Total Internal Reflection (TIR):**
  Geometrical ray bending and critical angle evaluation use the real part ``\text{Re}(\tilde{n}) = n``.
* **Optical Path Length (OPL):**
  The optical path length is computed as ``\text{OPL} = \int \text{Re}(\tilde{n}) \, \mathrm{d}s``, ensuring real-valued geometric path lengths.

---

## Refractive Index Types & Wrappers

Any callable object `f(λ::Real) -> Number` (returning a `Real` or `Complex` scalar) satisfies the [`RefractiveIndex`](@ref) interface. `BeamletOptics.jl` provides several built-in types:

```@docs; canonical=false
BeamletOptics.RefractiveIndex
```

### AbsorbingMedium

Use [`AbsorbingMedium`](@ref) to wrap a real base index with a known linear absorption coefficient $\alpha$ (in $[1/\text{m}]$):

```@docs; canonical=false
AbsorbingMedium
```

#### Example: Absorbing Neutral Density Filter / Lens
```julia
using BeamletOptics

# Absorbing medium with base index 1.5 and α = 100 1/m (1 1/cm)
med_abs = AbsorbingMedium(1.5, 100.0)

# Can also use dispersion functions:
med_disp = AbsorbingMedium(
    λ -> 1.51 + 0.005e-12 / λ^2,   # n(λ)
    λ -> 50.0 + 10.0e-6 / λ        # α(λ) in 1/m
)

# Create a 10 mm thick plano-planar absorbing block
slab = SphericalLens(Inf, Inf, 10e-3, 25.4e-3, med_abs)
```

### GainMedium

Use [`GainMedium`](@ref) to model active laser gain materials with small-signal gain coefficient $g$ (in $[1/\text{m}]$):

```@docs; canonical=false
GainMedium
```

#### Example: Nd:YAG Laser Crystal
```julia
using BeamletOptics

# Nd:YAG crystal with n = 1.82 and small-signal gain g = 50 1/m at 1064 nm
yag_gain = GainMedium(1.82, 50.0)

# 20 mm long AR-coated crystal
crystal = SphericalLens(Inf, Inf, 20e-3, 10e-3, yag_gain) |> with_coatings(
    front = SimpleARCoating(0.0),
    back  = SimpleARCoating(0.0)
)
```

### Sellmeier Dispersion

For transparent optical glasses specified by 6 Sellmeier coefficients:

```@docs; canonical=false
SellmeierEquation
```

### Tabulated / Interpolated Indices

For discrete tabulated refractive index data:

```@docs; canonical=false
DiscreteRefractiveIndex
```

### Direct Anonymous Functions

You can also directly pass anonymous functions returning complex numbers:

```julia
# Direct complex refractive index: n + i*κ
n_direct = λ -> complex(1.60, 5e-4)
lens = SphericalLens(50e-3, -50e-3, 5e-3, 25.4e-3, n_direct)
```

---

## Helper Functions

`BeamletOptics.jl` provides helper functions to query optical constants and compute attenuation factors:

```@docs; canonical=false
real_refractive_index
extinction_coefficient
absorption_coefficient
bulk_attenuation_factor
```

---

## Interaction with Optical Coatings

Bulk absorption and surface coatings operate independently and compose physically:
* **Bulk propagation:** Rays and beamlets attenuate as they travel through the medium according to $e^{-\alpha L}$.
* **Boundary physics:** When reaching a coated interface (e.g. [`ThinFilmCoating`](@ref) or [`SimpleARCoating`](@ref)), Fresnel reflection and transmission are evaluated at the boundary.
