# Tutorial draft: a pancake lens for a VR headset (roadmap #167, step 8)

> Draft for the docs, written against the **target API** of steps 5–8 (surface channels, coating
> models, coatings on faces, surface stacks). Names marked *(tentative)* do not exist yet. The draft
> fixes what the API must make easy; it is removed from the PR before merge and becomes a docs
> tutorial with `@example` blocks once the steps have landed.

## Why this showcase

A pancake lens folds the light path of a VR eyepiece three times inside a few millimetres. It only
works because **several thin layers act on one surface, in a fixed order, and the light passes
them several times**: a 50/50 mirror coated on a curved lens face, and a quarter-wave retarder
laminated under a reflective polarizer on the flat face. Every feature of steps 5–8 shows up:

| Feature | Where in the pancake lens |
|---|---|
| Channels and splitting (step 5) | the 50/50 mirror splits every ray, the polarizer sorts by polarization |
| Coating models (step 6) | 50/50 mirror, AR coating of the display cover |
| Coatings on faces (step 7) | mirror on the curved face, film stack on the flat face of the same lens |
| Surface stacks (step 8) | retarder + reflective polarizer, in physical order, passed twice |
| Non-sequential tracing | the solver finds the folded path by itself, no sequence is declared |
| Polarization | the fold is a polarization trick; a retarder error leaks a ghost image |

The result is a picture people remember: a ray bouncing back and forth in a lens 10 mm thick, and
a ghost image on the eye pupil that appears when the retarder is slightly off.

## 1. The folded path

```text
display ─▶ QWP₁ ─▶ 50/50 mirror (curved face) ─▶ glass ─▶ QWP₂ ─▶ reflective polarizer
                                                     ◀── reflected (wrong polarization)
           50/50 mirror, from inside ◀── glass ◀── QWP₂
                    reflected ──▶ glass ──▶ QWP₂ ──▶ reflective polarizer ──▶ eye (right polarization)
```

- The display emits linear polarization; QWP₁ makes it circular.
- Half of it passes the 50/50 mirror. QWP₂ turns it into the linear polarization that the
  reflective polarizer **reflects**.
- Back through QWP₂ it is circular again, the 50/50 mirror reflects half of it and flips its
  handedness, so after QWP₂ it is the linear polarization that the polarizer **transmits**.
- Ideal throughput: 0.5 × 0.5 = **25 %**, the textbook limit of a pancake lens. The folded path
  passes the curved mirror in reflection, which gives the optical power of the eyepiece.

## 2. Building the lens

```julia
using BeamletOptics, GLMakie

λ = 550e-9
n = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)

# Plano-convex lens: curved face towards the display, flat face towards the eye
lens = SphericalLens(Inf, -60e-3, 8e-3, 40e-3, n)

# 50/50 mirror on the curved face (step 7: faces by getters, a coating per face)
surface_model!(lens, back_face(lens), BeamsplitterCoating(0.5))          # (tentative names)

# Retarder under a reflective polarizer on the flat face, listed from the substrate outward
film = SurfaceStack(                                                      # (tentative)
    Retarder(π / 2, deg2rad(45)),       # quarter-wave plate, fast axis at 45°       (tentative)
    ReflectivePolarizer([1, 0, 0]),     # transmits x, reflects z                     (tentative)
)
surface_model!(lens, front_face(lens), film)
```

The order inside `SurfaceStack` is the physical order of the laminated films. A ray inside the
glass meets the retarder first, the polarizer second; the reflected channel of the polarizer passes
the retarder again on its way back. The user writes the films down as they are laminated; the
cascade of the surface model does the rest.

## 3. Display, eye and system

```julia
# Display: three field points, each a cone of polarized rays
display = [PointSource([x, -20e-3, 0], [0, 1, 0], deg2rad(15), λ; num_rays = 400,
                       polarization = [1, 0, 0]) for x in (-10e-3, 0, 10e-3)]  # (polarization kw tentative)
qwp1 = PolarizationFilter(30e-3, Retarder(π / 2, deg2rad(45)))               # (tentative)
translate3d!(qwp1, [0, -19e-3, 0])

translate3d!(lens, [0, 0, 0])
eye = Detector(8e-3)                     # the eye pupil, 15 mm eye relief
translate3d!(eye, [0, 23e-3, 0])

system = System([qwp1, lens, eye])
solve_system!(system, display)
```

## 4. What the user measures

1. **Throughput:** `optical_power(eye) / P_display` → 0.25 for ideal films. The test of the
   tutorial checks this against the analytic value.
2. **The folded path:** `render!(ax, system); render!(ax, display)` shows the rays reflected twice
   inside the lens. Showcase figure: a cut through the lens with the three passes coloured by
   pass number.
3. **The ghost:** the direct path (display → both faces → eye, without folds) is blocked by the
   polarizer only if QWP₁ and QWP₂ are exact. Scan the retardance error δ of QWP₂:

   ```julia
   for δ in deg2rad.(0:2:20)
       film = SurfaceStack(Retarder(π / 2 + δ, deg2rad(45)), ReflectivePolarizer([1, 0, 0]))
       surface_model!(lens, front_face(lens), film)
       solve_system!(system, display)
       # ghost: rays that reach the eye without reflection at the curved face
       ...
   end
   ```

   The leaked power grows as sin²(δ/2) of the input. The tutorial plots ghost-to-image ratio over δ
   and shows the ghost on the eye detector as an unfocused halo next to the sharp image: the reason
   why real pancake lenses need retarders with a tight tolerance and an achromatic design.
4. **Colour:** repeat at 450 nm and 650 nm. A zero-order QWP is a quarter wave at one wavelength
   only, so the ghost is coloured. This motivates dispersion in `Retarder` (an `AbstractRefractiveIndex`
   for the birefringence) as a later extension.

## 5. Requirements for steps 5–8 that this tutorial makes concrete

- A channel reflected by an outer layer passes the inner layers again (cascade, step 8).
- The 50/50 mirror on a body face acts from both sides; from inside the glass its reflection stays
  in the glass (enclosure of the reflected side, step 5).
- `AllChannels()` splitting with a power cutoff keeps the tree finite: each 50/50 pass halves the
  power, the cutoff stops the third and later bounces.
- A ray that ends at the reflective polarizer in the wrong polarization must not be lost when the
  power is checked: `optical_power` of the eye plus the power of rays that leave the system must
  equal the input, which is a good X1 test of the whole chain (energy conservation).
- A face getter for every catalog lens (`front_face`, `back_face`) and a setter of the face model.
