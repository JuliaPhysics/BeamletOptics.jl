# Polarizers

Polarizers act only on polarized light: `Beam(pos, dir, λ, E0)` (PolarizedRay beams) and
`AstigmaticGaussianBeamlet`s. Unpolarized `Ray`s pass an ideal filter unchanged.

| Constructor | Notes |
|-------------|-------|
| `PolarizationFilter(edge_length; cutoff_strength = eps())` | ideal, zero-thickness square in the x-z plane; **transmits x, blocks z** |
| `RoundPolarizationFilter(diameter; cutoff_strength = eps())` | round version |
| `RoundLinearPolarizer(diameter, front_thickness, back_thickness, n; cutoff_strength = eps())` | film between two glass plates: front glass at y ∈ [−front, 0], film at y = 0, back glass at y ∈ [0, back]. The uncoated glass surfaces add Fresnel losses. |
| `transmission_axis(p)` | current transmission axis (unit vector, global frame) |

- To set the angle, rotate the filter about the beam axis, e.g. `yrotate3d!(pf, θ)` for a +y beam.
  The sign of `transmission_axis` may flip; only its direction matters.
- A ray whose field falls below `cutoff_strength` is **terminated** at the filter. The beam then ends
  with the incident segment, so check `length(rays(beam))` before reading `last(rays(beam)).E0`.
- The model is a 3D Jones-matrix projection; out-of-plane tilts are handled by projection only.

See `templates/06_malus_polarizer.jl`.
