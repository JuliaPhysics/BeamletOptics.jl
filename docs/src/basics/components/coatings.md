# [Optical Coatings](@id coatings)

BeamletOptics.jl supports a flexible system for defining optical coatings on refractive and reflective components. This allows for the simulation of Anti-Reflection (AR) coatings, High-Reflectivity (HR) mirrors, beamsplitter coatings, and even spatially varying properties.

## Overview

The core of the system is the abstract type [`AbstractCoating`](@ref). All optical components that interact with rays (`Lens`, `Mirror`, `Prism`, etc.) now accept a `coating` parameter.

When a ray hits a surface, the `coating` determines the complex amplitude reflection ($r_s, r_p$) and transmission ($t_s, t_p$) coefficients based on the angle of incidence, wavelength, and refractive indices.

### Available Coating Types

*   [`Uncoated`](@ref): The default. Calculates standard Fresnel coefficients for a dielectric interface.
*   [`SimpleCoating`](@ref): Idealized coating defined by explicit Reflection ($R$) and Transmission ($T$) values (or functions).
*   [`MultilayerCoating`](@ref): Physically rigorous simulation of thin-film stacks using the Transfer Matrix Method (TMM).
*   [`SpatiallyVaryingCoating`](@ref): Applies different coatings to different parts of an object (e.g., separate coatings for the front, back, and side of a lens).

---

## Simple Coatings

For many simulations, you don't need to model the physical layer stack. You just want an "AR Coating" that reflects 0.1% or a Mirror that reflects 99%. Use [`SimpleCoating`](@ref) for this.

```julia
# 50/50 Beamsplitter coating (R=0.5, T=0.5)
bs_coating = SimpleCoating(0.5)

# Ideal Anti-Reflection coating (R=0.0)
ar_coating = SimpleCoating(0.0)

# High-Reflectivity mirror (R=0.99)
hr_coating = SimpleCoating(0.99)
```

You can also make them wavelength or angle dependent by passing functions:

```julia
# Reflectivity increases with angle
c = SimpleCoating(
    (λ, θ) -> 0.1 + 0.5 * sin(θ)^2,  # R(λ, θ)
    (λ, θ) -> 0.9 - 0.5 * sin(θ)^2   # T(λ, θ)
)
```

---

## Multilayer Coatings

For accurate physical modeling of interference effects, use [`MultilayerCoating`](@ref). This uses the **Transfer Matrix Method (TMM)** to calculate coefficients for arbitrary stacks of dielectric (or metallic) layers.

1.  Define your materials.
2.  Define the layer stack.
3.  Create the coating.

```julia
# Materials
SiO2 = RefractiveMaterial("SiO2", 1.45)
TiO2 = RefractiveMaterial("TiO2", 2.3)

# Design an antireflection coating for 550nm (Quarter-Wave stack)
λ0 = 550e-9
# Layers: Air | TiO2 (H) | SiO2 (L) | Glass
layers = [
    ThinFilmLayer(TiO2, λ0 / (4 * 2.3)),  # QWOT high index
    ThinFilmLayer(SiO2, λ0 / (4 * 1.45))  # QWOT low index
]

coating = MultilayerCoating(layers, λ0)
```

This coating will now automatically exhibit the correct angular and spectral dependence during ray tracing.

---

## Spatially Varying Coatings

Real optical components often have different coatings on different surfaces (e.g., AR on the front/back face, but absorbing or "ground" on the cylindrical edge).

The [`SpatiallyVaryingCoating`](@ref) type allows you to define a selection function that returns a specific coating based on the impact point `p` and surface normal `n`.

```julia
# Create a lens with AR coating on optical surfaces, but Absorbing on edges
front_ar = SimpleCoating(0.01) # 1% reflection
side_abs = SimpleCoating(0.0)  # Absorbing (R=0, T=0 if used on mirror/lens correctly)

# Note: BeamletOptics provides helper logic for this distribution automatically 
# when using the specific Lens constructors, but you can build one manually:

spatial_c = SpatiallyVaryingCoating((p, n) -> begin
    # Logic to decide which coating to return
    if abs(p[3]) < 0.05
        return front_ar
    else
        return side_abs
    end
end)
```

### Lens Specifics

The `Lens` constructor has been updated to accept separate arguments for `front_coating`, `back_coating`, and `side_coating` to make this easy:

```julia
lens = SphericalLens(...,
    front_coating = SimpleCoating(0.005), # High quality AR
    back_coating = SimpleCoating(0.005),
    side_coating = SimpleCoating(0.5)     # Scatter/Reflect on edge
)
```

Under the hood, this creates a `SpatiallyVaryingCoating` that uses exact distance checks against the lens geometry to ensure the correct coating is applied, even for complex shapes like meniscus lenses.

---

## API Reference

```@docs
AbstractCoating
Uncoated
SimpleCoating
MultilayerCoating
ThinFilmLayer
SpatiallyVaryingCoating
coating_coefficients
get_coating_at
```
