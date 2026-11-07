# Design: coatings on faces (roadmap #167)

> Design note for review, stacked on the coating models of #229. It is removed before the PR is merged.

## Decisions (maintainer, 2026-11-07)

1. **Faces are predicates:** `isonface(face, part, hit, point, normal)` in the local frame of the part.
   `FrontFace`, `BackFace`, `EdgeFace` (lens axis = local y), `FlatFace(normal)`, `AllFaces`. Only
   evaluated for parts with coated faces (`iscoated`, cached at the hit like `isbody`).
2. **The table hangs on the part:** `Lens` and `Prism` have a field `coatings::FaceModels`;
   `surface_model!(part, face, model)` sets, `surface_model(part, face)` reads, `nothing` removes. Other
   bodies opt in with `face_models(part)`.
3. **Masks are faces with an `AbstractShape`:** `MaskedFace(face, mask)`, inside via `isinside(shape,
   point)` (SDFs: `sdf ≤ 0`), the mask in the local frame of the part.
4. **System model:** `System(objects; surface_model = FresnelGhosts())`, also `StaticSystem`; default
   `Uncoated()`. Chain of `surface_model(t)`: thin layer → coated face of the body left → of the body
   entered → system.

## Boundary representations

A shape that knows its faces (BRep) sets the native face of each hit in `intersect3d` via
`face!(hit, face)`; the default `isonface` compares the face with the native face of the hit. Hence
BRep faces are coated directly with `surface_model!(part, brep_face, model)`, and geometric faces and
masks work on BReps too. `test/Components/TestFaces.jl` has a dummy shape with native faces as the
stand-in.

## Orientation

A coating lies on the outside of its body. `isfacing(t, ray)` tells a model whether the ray comes from
that side, relative to the normal of the surface that carries the model (thin layer, or the coated
body), not of the hit that the solver picked: at a cemented face the coating of one lens stays on its
own substrate.
