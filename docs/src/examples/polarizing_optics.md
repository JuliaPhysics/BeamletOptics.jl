# Polarizing Optics Examples

This page provides practical examples demonstrating how to use waveplates (Half-Wave Plates and Quarter-Wave Plates) and Polarizing Cube Beamsplitters in `BeamletOptics.jl` to manipulate and route polarized light.

---

## Circular Polarization Generation (Quarter-Wave Plate)

A Quarter-Wave Plate (QWP) introduces a relative phase retardance of $\pi/2$ (a quarter-wavelength) between the fast and slow axes of the waveplate. When linearly polarized light passes through a QWP at an angle of 45° relative to the fast axis, the output light becomes circularly polarized.

The following example demonstrates how to set up a flat, circular QWP to generate circularly polarized light:

```@example polarizing_showcase
using BeamletOptics
using LinearAlgebra
const BMO = BeamletOptics

# Define QWP (diameter = 10 mm, retardance = π/2) with fast axis at 0° (along x-axis)
qwp = QuarterWavePlate(10e-3; fast_axis_angle = 0.0)
system = System([qwp])

# Input beam linearly polarized at 45° in the transverse plane (x-z plane)
# propagation along y-axis
ray_in = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, [1.0, 0.0, 1.0] / sqrt(2))
beam = Beam(ray_in)

# Trace ray through system
solve_system!(system, beam; retrace = false)

# Analyze exit polarization
ray_out = rays(beam)[2]
E_out = BMO.polarization(ray_out)

println("Input polarization:  ", BMO.polarization(ray_in))
println("Output polarization: ", E_out)
println("Magnitude X-component: ", abs(E_out[1]))
println("Magnitude Z-component: ", abs(E_out[3]))
```

As shown above, the $x$- and $z$-components of the output electric field have equal magnitude ($\approx 0.707$) and a phase difference of 90° ($e^{i \pi/2} = i$), which corresponds to circularly polarized light.

---

## Linear Polarization Rotation (Half-Wave Plate)

A Half-Wave Plate (HWP) introduces a relative phase retardance of $\pi$ (a half-wavelength) between the fast and slow axes. If linearly polarized light is incident on a HWP with its polarization plane at an angle $\theta$ relative to the fast axis, the plane of polarization is rotated by $2\theta$.

In this example, we rotate horizontal linear polarization ([1, 0, 0]) by 45° by setting the HWP's fast axis at 22.5°:

```@example polarizing_showcase
# Create flat HWP with fast axis oriented at 22.5° (0.3927 radians)
hwp = HalfWavePlate(10e-3; fast_axis_angle = deg2rad(22.5))
system_hwp = System([hwp])

# Input horizontal polarization along x-axis
ray_in_hwp = PolarizedRay([0.0, -10e-3, 0.0], [0.0, 1.0, 0.0], 1000e-9, [1.0, 0.0, 0.0])
beam_hwp = Beam(ray_in_hwp)

# Trace ray through system
solve_system!(system_hwp, beam_hwp; retrace = false)

# Output polarization
ray_out_hwp = rays(beam_hwp)[2]
E_out_hwp = BMO.polarization(ray_out_hwp)

println("Input polarization:  ", BMO.polarization(ray_in_hwp))
println("Output polarization: ", E_out_hwp)
```

The output polarization vector is approximately `[0.0 - 0.707im, 0.0, 0.0 - 0.707im]`, which corresponds to a linear polarization state oriented at 45° in the transverse plane with an overall global phase shift of -90° (represented by the common factor of `-im`).

---

## Polarization Routing (Polarizing Beamsplitter Cube)

A Polarizing Beamsplitter Cube (PBS) splits an incoming beam into two polarized components. Light polarized perpendicular to the plane of incidence (s-polarization) is reflected at a 90° angle, while light polarized parallel to the plane of incidence (p-polarization) is transmitted straight through the cube.

This example splits an input beam containing equal s- and p-components:

```@example polarizing_showcase
# Define PBS with leg length 25 mm and refractive index n = 1.5
pbs = PolarizingCubeBeamsplitter(25e-3, n -> 1.5)
# Shift it forward so the center lies at y = 50 mm
translate3d!(pbs, [0.0, 50e-3, 0.0])
system_pbs = System([pbs])

# Input beam with equal s- and p-polarized components (linear at 45°)
ray_in_pbs = PolarizedRay([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1000e-9, [1.0, 0.0, 1.0] / sqrt(2))
beam_pbs = Beam(ray_in_pbs)

# Trace the beam (depth_max = 2 to allow multiple reflections/refractions through the cube boundaries)
solve_system!(system_pbs, beam_pbs; depth_max = 2, retrace = false)

# Examine child beams
# children[1] is the transmitted beam
# children[2] is the reflected beam
transmitted = beam_pbs.children[1]
reflected = beam_pbs.children[2]

println("Number of child beams: ", length(beam_pbs.children))
println("Transmitted beam direction: ", direction(last(rays(transmitted))))
println("Transmitted polarization:     ", BMO.polarization(last(rays(transmitted))))
println("Reflected beam direction:   ", direction(last(rays(reflected))))
println("Reflected polarization:       ", BMO.polarization(last(rays(reflected))))
```

As demonstrated:
* The transmitted beam propagates straight along `[0, 1, 0]` and contains the p-polarized component scaled to $\approx [0.679, 0.0, 0.0]$ due to Fresnel reflection losses at the entrance and exit surfaces of the glass cube.
* The reflected beam is redirected by 90° along `[-1, 0, 0]` and contains the s-polarized component $\approx [0.0, 0.0, -0.679]$ (where the negative sign arises from local coordinate conventions).
