# Polarization optics

Polarization sensitive components use the `PolarizedRay` type and the Jones formalism
implemented in this package. The following constructors create simple waveplates
and polarizing beamsplitters.

```@docs
Waveplate
RoundWaveplate
QuarterWaveplate
HalfWaveplate
RectangularPolarizingPlateBeamsplitter
RoundPolarizingPlateBeamsplitter
PolarizingThinBeamsplitter
RoundPolarizingThinBeamsplitter
PolarizingCubeBeamsplitter
```

## Example

The snippet below demonstrates how a half-wave plate in combination with a
polarizing cube beamsplitter can be used to rotate linear polarization.

```@example polarization_example
using BeamletOptics
mm = 1e-3
wp = HalfWaveplate(12mm)
cbs = PolarizingCubeBeamsplitter(20mm, λ->1.5)
translate3d!(wp, [0, 30mm, 0])
system = StaticSystem([cbs, wp])
ray = PolarizedRay([0,0,0], [0,1,0], 532e-9, [1,0,0])
beam = Beam(ray)
solve_system!(system, beam)
polarization(first(rays(beam.children[1])))
```
