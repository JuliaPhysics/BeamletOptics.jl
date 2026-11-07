# Coatings

A coating is a surface model that changes reflection and transmission at a boundary. Coat a **face of a
body** (`Lens`, `Prism`, so also plates, cube halves and the lenses of a doublet) with
`surface_model!(part, face, model)`; `surface_model(part, face)` reads it, `surface_model!(part, face,
nothing)` removes it. The coating moves with the part and acts between the media on both sides, also at a
cemented face.

| Face | Meaning |
|------|---------|
| `FrontFace()`, `BackFace()`, `EdgeFace()` | lens faces: outward normal against, along, perpendicular to the lens axis (local +y) |
| `FlatFace(normal)` | flat face of a prism or plate with this outward normal in the **local** frame of the part |
| `AllFaces()` | every face |
| `MaskedFace(face, mask)` | the part of `face` inside the shape `mask` (local frame), e.g. `BMO.CylinderSDF(r, h)` (along local y), `BMO.RingSDF`, `BMO.BoxSDF` |

If several coated faces apply to a hit, the one coated first acts. A free-standing film (pellicle, film
between two prisms) is a `CoatingLayer(shape, model)` with its own 2D shape; its front is the side its
normal points to (`RectangularFlatMesh`: +y). `System(objects; surface_model = FresnelGhosts())` (also
`StaticSystem`) sets the model of every uncoated boundary.

```julia
using BeamletOptics
const BMO = BeamletOptics

λ0 = 550e-9
lens = SphericalLens(50e-3, -50e-3, 5e-3, 25e-3, ConstantRefractiveIndex(1.52))
ar = ThinFilmCoating(QuarterWaveLayer(ConstantRefractiveIndex(1.38), λ0); splitting = BMO.Dominant())
surface_model!(lens, FrontFace(), ar)
surface_model!(lens, BackFace(), ar)
surface_model!(lens, MaskedFace(FrontFace(), BMO.CylinderSDF(3e-3, 10e-3)), MirrorCoating())  # coated first wins: AR
translate3d!(lens, [0, 50e-3, 0])

# A free-standing longpass dichroic at 45°: reflects below 800 nm, transmits above
dichroic = CoatingLayer(BMO.RectangularFlatMesh(25e-3, 25e-3),
    SpectralCoating((λ, θ) -> λ < 800e-9 ? 0.98 : 0.02))
zrotate3d!(dichroic, deg2rad(45))
translate3d!(dichroic, [0, 100e-3, 0])

beam = Beam([0, 0, 0], [0, 1, 0], 532e-9)
solve_system!(System([lens, dichroic]), beam)
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
- Every split doubles the beams: bound the tree with `depth_max` of `solve_system!` or the preference
  `power_cutoff` (`get_power_cutoff()`, fraction of the root power; set via
  `BeamletOptics.Config.set_preference!("power_cutoff", 1e-4)` and restart Julia).
- The lossless coatings reflect with an extra π when the beam comes from the front of the coating, i.e.
  from outside of the coated body or from the side the layer normal points to (Stokes relations);
  a `ThinFilmCoating` takes its phases from the design.
- Beams without polarization (`GaussianBeamlet`) carry the power fractions of unpolarized light.
