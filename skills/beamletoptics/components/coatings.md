# Coatings

A coating is a surface model that changes reflection and transmission at a boundary. Place it as a
`CoatingLayer(shape, model)`: a thin layer flush on the face of a substrate (the solver combines layer and
substrate at the same boundary) or free-standing. The **normal of the layer shape points away from the
substrate**; `RectangularFlatMesh` has its normal along +y, so rotate it with `zrotate3d!(layer, π)` for the
front face of an element that faces −y.

```julia
using BeamletOptics
const BMO = BeamletOptics

# A free-standing longpass dichroic at 45°: reflects below 800 nm, transmits above
dichroic = CoatingLayer(BMO.RectangularFlatMesh(25e-3, 25e-3),
    SpectralCoating((λ, θ) -> λ < 800e-9 ? 0.98 : 0.02))
zrotate3d!(dichroic, deg2rad(45))
translate3d!(dichroic, [0, 50e-3, 0])

# A single-layer AR coating; flush on a face, the layer must lie within 1 nm of it
ar = ThinFilmCoating(QuarterWaveLayer(ConstantRefractiveIndex(1.38), 550e-9); splitting = BMO.Dominant())

beam = Beam([0, 0, 0], [0, 1, 0], 532e-9)
solve_system!(System([dichroic]), beam)
```

| Model | Use |
|-------|-----|
| `BeamsplitterCoating(r, t)` | lossless split with field amplitudes `r`, `t` (`r² + t² = 1`), carried by `ThinBeamsplitter` |
| `MirrorCoating()` | ideal mirror, `R = 1` |
| `AntireflectionCoating()` | ideal AR, `R = 0`, no Fresnel loss |
| `SpectralCoating(R)` / `SpectralCoating(Rs, Rp)` | power reflectance `R(λ, θ)` (λ in m, θ in rad), s and p separately; `SpectralCoating(λs, Rs)` interpolates a table linearly, constant beyond it |
| `ThinFilmCoating(layers...)` | coherent thin-film stack, layers `ThinFilmLayer(n, d)` or `QuarterWaveLayer(n, λ0)` **from the substrate outward**; absorbing layers via `AbsorbingRefractiveIndex` |
| `FresnelGhosts()` | the bare boundary with both Fresnel parts traced (ghost analysis) |

- Splitting coatings (`BeamsplitterCoating`, `SpectralCoating`, `ThinFilmCoating`) trace the reflected and
  the transmitted light as `children` `[transmitted, reflected]` by default; `splitting = BMO.Dominant()`
  (keyword of `SpectralCoating` and `ThinFilmCoating`) traces only the stronger part, e.g. for AR coatings.
- Every split doubles the beams: bound the tree with the keywords `power_cutoff` (fraction of the root
  power, e.g. `solve_system!(system, beam; power_cutoff = 1e-4)`) or `depth_max` of `solve_system!`.
  Beams of `Ray`s carry no power, only `depth_max` bounds them.
- The lossless coatings reflect with an extra π when the beam faces the layer normal (Stokes relations);
  a `ThinFilmCoating` takes its phases from the design.
- Beams without polarization (`GaussianBeamlet`) carry the power fractions of unpolarized light.
