```@setup detectors
detector_showcase_dir = joinpath(@__DIR__, "..", "..", "assets", "detector_assets")

Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "spotdetector_showcase.jl"), use_placeholder=false)
Main.DocUtils.conditional_include(joinpath(detector_showcase_dir, "photodetector_showcase.jl"), use_placeholder=false)
Main.DocUtils.conditional_include(joinpath(@__DIR__, "..", "..", "assets", "examples", "psfdetector_showcase.jl"), use_placeholder=false)
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

### Point spread function estimation

If a [`Detector`](@ref) is placed in the focal plane of an imaging system, the coherent addition of the ray-attached plane waves yields an estimate of its point spread function (PSF). The call is the same as above, only the detector position changes. A singlet imaging a collimated 15 mm beam yields the expected Airy pattern:

![Airy disc PSF](psf_airy_showcase.png)

For [`PolarizedRay`](@ref)s the field is added as 3D vectors, so [`electric_field`](@ref) returns a matrix of complex vectors and each component of the focal field can be analyzed separately. This becomes relevant at high NA, where the field vectors tilt towards the optical axis.

!!! warning "Experimental feature"
    The PSF estimation does not use pupils (yet) but merely superimposes the ray-attached plane waves. It gives qualitatively sound results, but requires good sampling of the problem to be quantitatively meaningful. No Strehl ratio is calculated.

!!! tip "PSF examples"
    The [Point spread functions](@ref) example page covers the Airy disc, an aberrated asphere and the vectorial focus of a parabolic mirror at NA 0.88, including the caveats on collimated sources and ray amplitudes.
