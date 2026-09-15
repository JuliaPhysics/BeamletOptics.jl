```@setup detectors
detector_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "detector_assets")

Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "spotdetector_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "photodetector_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "psfdetector_showcase.jl"))
Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "polarized_psf_showcase.jl"))
```

# Detectors

In practice, photodetectors allow the conversion of electromagnetic radiation into electric signals. BMO provides the [`Detector`](@ref) as an element to capture and evaluate optical data during and after running a simulation, respectively. The `Detector` is designed to accumulate e.g. field or ray data, enabling analysis of intensity distributions, interference patterns, and other beam properties. Currently, post-processing capabilities are limited to the functionality as described in the [Spot diagrams](@ref) and [Field distributions](@ref) sections.

In general, detector-like elements are supposed to fall under the [`BeamletOptics.AbstractDetector`](@ref) type, which defines a interface for detector implementations.

!!! warning "Resetting detectors"
    In general, the data stored in a `Detector` is not automatically reset between calls of [`solve_system!`](@ref). This task is placed within the responsibility of the user. A detector reset can be performed with the [`empty!`](@ref) function.

## Detector type

A `Detector` can be easily spawned by initializing e.g. `pd = Detector(5mm)` which will create a 5x5 mm² detection screen.

```@docs; canonical=false
Detector(::Real, ::Bool)
Detector
```

After solving a system containing a `Detector`, the methods listed below can be used in order to analyze the stored data. If no data is obtained during the tracing procedure, an error message will be stored.

## Spot diagrams

The `spot_diagram` method provides a straight forward way to generate spot diagrams, which are commonly used to perform initial assessments of the optical performance of an imaging setup.

```@docs; canonical=false
spot_diagram
```

Below an optical system consisting of a collection of collimated [`Beam`](@ref)s passing through a [`ThinLens`](@ref) is shown. A [`Detector`](@ref) is positioned at the approximate focal plane to capture the resulting spot diagram.

![Thin lens setup](spot_diagram_system.png)

The beam bundle used to generate the spot diagram was created via the [`CollimatedSource`](@ref) constructor. The resulting spot diagram of the lens shown above is visualized below.

![Spot diagram showcase](spot_diagram_showcase.png)

## Field distributions

Alternatively, the detector data can also be used to reconstruct electric field distributions of incoming beams on its surface using coherent addition. Depending on the beam type, either plane wave or Gaussian beam models are used. As a user, this data can be accessed using the [`electric_field`](@ref) interface.

```@docs; canonical=false
electric_field(::Detector)
```

For convienience, the [`intensity`](@ref) function returns flux values directly. The [`optical_power`](@ref) method can be used in order to obtain the total optical power on the detector surface.

```@docs; canonical=false
intensity(::Detector)
```

### Gaussian beamlet interference

One of the use cases of the [`Detector`](@ref) is to analyse interference patterns. Below a rendered example of a detector model ([FDS010](https://www.thorlabs.com/thorproduct.cfm?partnumber=FDS010)) can be seen. The detector active area is marked in blue (1x1 mm²). 

![Photodetector showcase](pd_showcase.png)

The figure below demonstrates an example intensity distribution captured by the detector pictured above, showing radial fringes due to a mismatch of the radii of curvature of the interfering [`GaussianBeamlet`](@ref)s. In addition, the beam waists have been visualized.

!!! tip "Interferometer tutorial"
    Refer to the [Michelson interferometer](@ref) section for a detailed tutorial on how to use the [`Detector`](@ref).

![Interference fringes showcase](fringes_showcase.png)

### Point spread function calculation

!!! warning "Experimental feature"
    The point spread function estimation is a highly experimental feature. It does not use
    pupils (yet) but merely uses superposition of the ray-attached plane-waves. While this
    gives qualitatively sound results, it requires good sampling of the problem to obtain
    quantitatively good results. Currently no Strehl-ratio is calculated due to that.

The package offers a simple method to estimate the point spread function of a system. It is 
currently limited and requires careful assessment by the user, if the results are to be trusted.

To analyze the PSF of a imaging system a [`Detector`](@ref) is added to the system at the plane
and orientation where the PSF is requested. This is the same approach as for the other detector types. 
The intensity map together with the coordinate system of the detector can be retrieved after solving
the system by calling the [`intensity`](@ref) function.

!!! warning
    When dealing with a collimated source as the input to your optical system, where you want to calculate the PSF, **DO NOT** use the [`CollimatedSource`](@ref) beam group directly but instead use the [`UniformDiscSource`](@ref) constructor. This function returns a `CollimatedSource` with an equal-area sampling, which correctly weights the outer beams in relation to the inner beams. Otherwise the results might be wrong.

#### Airy-Disc Example

This is a classic example where a collimated circular beam is imaged onto a point by a singlet lens.
Due to the finite size of the aperture stop (in this case given by the 15 mm size of the beam), the diffraction
limited intensity pattern is given by the Airy-disc:

```math
I(r)=I_0\!\left[\frac{2J_1\!\bigl(\pi D r/(\lambda f)\bigr)}{\pi D r/(\lambda f)}\right]^2
```

With ``r`` the radius from the origin, ``I_0`` the maximum intensity, ``J_1`` the Bessel function of the first kind of order one, ``D`` the aperture width, ``\lambda`` the wavelength and the focal length ``f``.

```julia
# example parameters
l = 1e-3
R1 = 100e-3
R2 = Inf
d = 25.4e-3
n = 1.5
λ = 1e-6

# generate uniform source, lens and detector
cs = UniformDiscSource([0, -10mm, 0], [0, 1, 0], 15e-3, λ)
lens = SphericalLens(R1, R2, l, d, x -> n)
detector = Detector(10e-3)

# shift detector into focus
translate3d!(detector, [0, 200mm + 0.13mm, 0])

# build system
sys = System([lens, detector])

solve_system!(sys, cs)

# retrieve intensity
x, z, I_num = intensity(detector; n=500, crop_factor=10)
```

Visualizing the result yields the expected Airy-disk pattern.

![Airy disc PSF](psf_airy_showcase.png)

#### Coma and Astigmatism Example

In this example, an aspheric lens images the collimated source onto a point but is tilted around the x-axis by 0.5 degrees.
This results in aberrations distorting the stigmatic imaging and leading to coma and astigmatism.

```julia
k = -0.675
d = 75.0e-3
l = 15e-3
radius = 76.68e-3
A = [0*(1e3)^1, 2.7709219e-8*(1e3)^3, 6.418186e-13*(1e3)^5, -1.5724014e-17*(1e3)^7, -2.7768768e-21*(1e3)^9, -2.590162e-25*(1e3)^11]
AL75150 = Lens(
    EvenAsphericalSurface(radius, d, k, A),
    l,
    n -> 1.5006520430
)

xrotate3d!(AL75150, deg2rad(-0.5))

detector = Detector(15e-3)

translate3d!(detector, [0, 158.1779e-3, 0.0])
system = System([AL75150, detector])

ps = UniformDiscSource([0, -0.1, 0], [0,1,0], 0.8*d, 1550e-9)

solve_system!(system, ps)
```

![Tilted asphere PSF](psf_tilted_showcase.png)

#### Vectorial focusing (high NA)

When a system is traced with [`PolarizedRay`](@ref)s, the detector adds the field vectors `E0` of all rays coherently as 3D vectors instead of scalar phasors. For small numerical apertures (NA) this reproduces the scalar result, since all field vectors are nearly parallel. At high NA, however, the rays converge from steep angles and their field vectors tilt towards the optical axis. As described by Richards and Wolf [Richards1959](@cite), this has two visible consequences for a linearly polarized input:

- the PSF is stretched along the direction of the input polarization, and
- a longitudinal field component (along the optical axis) appears, which vanishes on axis but forms two lobes along the polarization direction.

BMO does not provide an ideal aplanatic lens, so this example constructs the converging ray fan of an aplanatic system with NA = 0.9 directly on a reference sphere around the focus. The optical axis is ``+y``, the input polarization is ``x``, and each ray carries the Richards–Wolf field and apodization ``\sqrt{\cos\alpha}``, weighted by the ring area element ``\sin\alpha``.

!!! warning "Internal API"
    The detector hits below are created manually via non-exported functions (`BeamletOptics.PolarizedRayHit`, `BeamletOptics.intersection!`) to emulate an ideal focusing element. In a regular simulation hits are generated by [`solve_system!`](@ref).

```julia
using BeamletOptics, LinearAlgebra
const BMO = BeamletOptics

λ = 1e-6

function focus_hit(dir, E0; focus = zeros(3))
    ray = PolarizedRay(focus .- λ .* dir, dir, λ, E0)
    BMO.intersection!(ray, BMO.Intersection(λ, BMO.Point3(-dir)))
    return BMO.PolarizedRayHit(ray, BMO.optical_path_length(ray))
end

function aplanatic_focus(NA; N = 60, pol = [1.0, 0, 0])
    pd = Detector(1.0)
    for i in 1:N, j in 1:(4N)
        α = asin(NA) * (i - 0.5) / N
        φ = 2π * (j - 0.5) / (4N)
        dir = [sin(α) * cos(φ), cos(α), sin(α) * sin(φ)]
        e_r = [cos(φ), 0, sin(φ)]
        e_φ = [-sin(φ), 0, cos(φ)]
        e_ρ = [cos(α) * cos(φ), -sin(α), cos(α) * sin(φ)]
        E0 = sqrt(cos(α)) * sin(α) .* (dot(pol, e_r) .* e_ρ .+ dot(pol, e_φ) .* e_φ)
        push!(pd, focus_hit(dir, E0))
    end
    return pd
end

NA = 0.9
R = 1.5λ / NA
pd = aplanatic_focus(NA)
xs, zs, E = electric_field(pd; n = 201, x_min = -R, x_max = R, z_min = -R, z_max = R)

I  = intensity.(E)                # total intensity
Ey = map(e -> abs2(e[2]), E)      # longitudinal component
```

The field `E` is a matrix of complex 3D vectors, so individual components can be evaluated separately. The cuts through the focus show a full width at half maximum that is about 1.3 times larger along ``x`` than along ``z``, and the longitudinal component reaches roughly 15 % of the peak of the transverse component.

![Vectorial PSF at NA 0.9](psf_vector_showcase.png)

!!! note "Unpolarized light"
    Orthogonal field vectors do not interfere. To obtain the PSF of unpolarized light, trace the system once with `pol = [1, 0, 0]` and once with `pol = [0, 0, 1]` and add the resulting intensities.