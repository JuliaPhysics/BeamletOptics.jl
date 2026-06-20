# Coatings

In `BeamletOptics.jl`, coatings are physical models attached to optical boundaries. They modify the complex amplitude and polarization state of propagated rays and beamlets (such as `PolarizedRay` or `AstigmaticGaussianBeamlet`) according to their physical structure.

The coatings system is designed around a clean separation between **geometry** (the shape of a lens, mirror, or flat interface) and **physics** (how light interacts with the boundary).

---

## Coating Behaviors

Every coating has an associated `CoatingBehavior` trait. The core ray-tracing engine uses this trait to determine the control flow of the ray propagation (e.g., whether a ray transmits, reflects, or splits):

* **`Transmissive`**: Light is refracted through the interface into the adjacent medium. (Default for `Uncoated`, `SimpleARCoating`).
* **`Reflective`**: Light is reflected back into the incident medium. (Default for `SimpleHRCoating`).
* **`Splitting`**: The incident ray is explicitly split into both a transmitted ray and a reflected ray. (Default for `SimpleBeamsplitterCoating`, or when modeling beamsplitter interfaces).
* **`Absorptive`**: Light is immediately absorbed (terminated) at the interface, generating no child rays or beamlets.

### Custom / Dynamic Coating Behaviors

By default, the core engine retrieves the behavior trait of a coating statically by calling `coating_behavior(coating)`. However, developers can implement custom and dynamic behaviors by overloading the two-argument version:

```julia
BeamletOptics.coating_behavior(coating::MyCustomCoating, ray::AbstractRay)
```

This allows a coating to dynamically change its behavior (e.g., switching from `Transmissive()` to `Reflective()` or `Splitting()`) on-the-fly based on properties of the incident ray—such as its angle of incidence, wavelength, or polarization state.

For example, a custom wavelength-selective coating could be implemented as:
```julia
function BeamletOptics.coating_behavior(coating::MyFilter, ray::AbstractRay)
    if wavelength(ray) < 600e-9
        return Reflective()
    else
        return Transmissive()
    end
end
```

---


## Available Coating Models

`BeamletOptics.jl` provides several built-in coating models, ranging from idealized parameters to exact physical wave-propagation models:

### Uncoated
Represents a bare interface between two media. The Fresnel transmission and reflection coefficients are computed dynamically based on the angle of incidence and the refractive indices of the adjacent media.
```@docs; canonical=false
BeamletOptics.Uncoated
```

### Simple AR / HR Coating
Idealized anti-reflective (AR) or highly-reflective (HR) coatings defined by a constant power reflectance $R$.
```@docs; canonical=false
BeamletOptics.SimpleARCoating
BeamletOptics.SimpleHRCoating
```

### Simple Beamsplitter Coating
An idealized splitting coating defined by constant complex amplitude coefficients for reflection ($r_s, r_p$) and transmission ($t_s, t_p$).
```@docs; canonical=false
BeamletOptics.SimpleBeamsplitterCoating
```

### Jones Coating
A generalized polarizing and phase-shifting boundary model defined by custom Jones matrices for transmission and reflection. The matrices can be static or functions of wavelength $\lambda$.
```@docs; canonical=false
BeamletOptics.JonesCoating
```

### Thin-Film Interference Coating
A physical multi-layer thin-film stack defined by a vector of refractive indices $n_j(\lambda)$ (which can be complex to model absorbing/metallic media) and physical thicknesses $d_j$. It calculates exact amplitude and phase coefficients using the **characteristic transfer matrix method**.
```@docs; canonical=false
BeamletOptics.ThinFilmCoating
```

---

## Applying Coatings to Components

You can attach coatings to standard optical components (like lenses or mirrors) using the fluent `with_coatings` API or by directly constructing `CoatedRefractive` or `CoatedMirror` wrappers.

```@docs; canonical=false
BeamletOptics.with_coatings
BeamletOptics.CoatedRefractive
BeamletOptics.CoatedMirror
```

### Face-Selective Coating
When coating a lens or prism, you can target specific faces by passing a `Pair` containing a filter pattern and the coating model.

The filter pattern can be:
* A `Symbol` representing a named face of the shape (e.g., `:front`, `:back`, `:side` for rotationally symmetric lenses, or `:hypotenuse`, `:leg1`, `:leg2` for a `RightAnglePrism`).
* An orientation vector (`AbstractVector`), matching faces whose local normal vector has a positive dot product with the orientation vector.
* A spatial predicate function `(local_pos, local_normal) -> Bool`.

#### Examples
```julia
# Coat only the front of a lens
coated_lens = lens |> with_coatings(front = SimpleARCoating())

# Coat the hypotenuse of a prism with HR, and Leg 1 with AR
coated_prism = with_coatings(prism, :hypotenuse => SimpleHRCoating(), :leg1 => SimpleARCoating())

# Spatial predicate: coat only the left half (x < 0) of a lens surface
left_half_ar = with_coatings(lens, ((p, n) -> p[1] < 0) => SimpleARCoating())
```

### Direct Construction

Instead of using the `with_coatings` API, you can directly construct the wrapper types `CoatedRefractive` (which also has the alias `CoatedLens`) and `CoatedMirror` using keyword arguments or explicit coating tuples.

#### Using Keyword Arguments
For convenience, you can pass `front` and `back` keyword arguments to attach coatings directly to the front/back surfaces:
```julia
# Direct construction of a coated lens
coated_lens = CoatedLens(lens; front = SimpleARCoating(), back = SimpleARCoating())

# Direct construction of a coated mirror
coated_mirror = CoatedMirror(mirror; front = SimpleHRCoating())
```

#### Using Explicit Coating Tuples
You can also construct the wrappers by passing a tuple of filter-to-coating-model pairs. This is particularly useful for complex multi-face shapes (like prisms) or when programmatically applying coatings:
```julia
# Attach coatings to specific faces of a prism via mapping tuples
coated_prism = CoatedRefractive(prism, (
    :hypotenuse => SimpleHRCoating(),
    :leg1 => SimpleARCoating()
))
```

---

## Standalone Coatings

If you want to represent a thin, free-floating coated interface (such as a flat filter or beam-splitting window) without wrapping a bulk refractive optic, you can use the standalone `Coating` component. This wraps a flat 2D shape and a coating model.

```@docs; canonical=false
BeamletOptics.Coating
```
