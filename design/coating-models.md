# Design: coating models (roadmap #167, step 6)

> Design note for review, stacked on the surface channels of #228. It is removed before the PR is
> merged; what survives goes into the docstrings, a Basics page on coatings and "Surface channels" in
> `docs/src/api/core.md`.

## Goal

Real coatings on real optics, as surface models on the channel contract of step 5: anti-reflection
and mirror coatings, beamsplitters and dichroics that depend on wavelength, angle and polarization,
and a coherent thin-film model for a coating whose layers are known. Every model only implements
`channels(model, t, ray, cs)` (and `splitting_of`); none has code per beam type.

Until faces can carry models (step 7), a coating acts as a **thin layer**: an object with a shape and
a model, which the coincident-hit resolution of #177/#178 already places on a substrate. Step 7 then
assigns the same model to a face; the models of this step do not change.

## Models

| Model | Physics | Channels | `splitting_of` |
|---|---|---|---|
| `BeamsplitterCoating(R)` (exists) | lossless, fixed `R` | `r`, `t`, π when facing the normal | `AllChannels` |
| `AntireflectionCoating()` | ideal AR: no reflection | transmitted with `|t| = 1` power-normalized, reflected dark | `Dominant` |
| `MirrorCoating()` | ideal HR: `R = 1` | reflected with π when facing the normal, transmitted dark | `Dominant` |
| `SpectralCoating(Rs, Rp)` | `Rs(λ, θ)`, `Rp(λ, θ)` given as functions or tables, lossless or with absorption `A` | s and p separately: `SPBasis(±√Rs, 0, 0, ∓√Rp)`, `SPBasis(√Ts, 0, 0, √Tp)` | `AllChannels` |
| `ThinFilmStack(layers, n_ref)` | coherent layers by the characteristic matrix (Born & Wolf §1.6, Macleod ch. 2), complex `rs, rp, ts, tp` for both incidence sides | Fresnel-like `SPBasis` from the stack | `AllChannels` |
| `FresnelGhosts()` | the bare boundary, both channels traced (decision 2 of #228) | the channels of the bare boundary, scalar `a` from the unpolarized power | `AllChannels` |

- **`SpectralCoating`** covers dichroics, edge filters and datasheet curves; the OpenRaman longpass
  mirror of the tutorials becomes one line instead of an own `interact3d`.
- **`ThinFilmStack`** computes everything from the layers: AR, HR, beamsplitters and filters follow
  from their designs, including the phase, so the convention of `conventions.md` is a consequence,
  not a setting. Layers are `(n::AbstractRefractiveIndex, d)` from the substrate outward, the same
  physical order as the stacks of step 8.
- **Several thin layers at one boundary** (today only one of them acts, see #178) combine into one
  model: a `ThinFilmStack` if they are thin films, otherwise the cascade of step 8.

## A thin layer that carries any model

`ThinBeamsplitter` is a thin layer with a `BeamsplitterCoating`. A generic thin layer carries any
model, so a user can put a coating onto a substrate today:

```julia
layer = CoatingLayer(RectangularFlatMesh(25e-3, 25e-3), ThinFilmStack(...))   # (tentative)
```

`ThinBeamsplitter(shape, R)` stays as the catalog name for the common case.

## Tests (X1)

- Quarter-wave AR on BK7 at its design wavelength: `R` against `(n0 ns - n1²)² / (n0 ns + n1²)²`.
- Quarter-wave HR stack of `N` pairs: `R` against the analytic formula, rising with `N`.
- Lossless stacks: `R + T = 1` for s and p at oblique incidence; reciprocity of `T` for both sides.
- `ThinFilmStack` without layers equals the bare Fresnel boundary.
- A polarized Mach-Zehnder with `ThinFilmStack` beamsplitters conserves power (the phase follows from
  the stack).
- `SpectralCoating` from a table reproduces the OpenRaman longpass result.
- `FresnelGhosts` on a plate: the first ghost carries `R²` of the root power.

## Open questions for the maintainers

1. **Absorbing layers** need complex indices, but the contract of `AbstractRefractiveIndex` returns a
   `Real` (#213). Extend it with an `extinction_coefficient(n, λ)` (default 0, step 11 of #167) or
   allow a complex return value for coating layers only?
2. **Names:** `AntireflectionCoating`, `MirrorCoating`, `SpectralCoating`, `ThinFilmStack`,
   `FresnelGhosts`, `CoatingLayer`. Catalog style `[Shape][Variant][Family]` suggests the family
   `Coating` as suffix; `ThinFilmStack` and `FresnelGhosts` are models, not catalog parts.
3. **`SpectralCoating` input:** functions `R(λ, θ)` only, or also tables with interpolation? Tables
   need an interpolation dependency or a small linear interpolation in BMO.
4. **`FresnelGhosts` as default model of a system** (`surface_model(system)`): in this step or with
   the faces in step 7?
