```@setup plots
using GLMakie, BeamletOptics

const BMO = BeamletOptics

tutorial_dir = joinpath(@__DIR__, "..", "assets", "ms_assets")

include(joinpath(tutorial_dir, "miniscope_showcase.jl"))
include(joinpath(tutorial_dir, "miniscope_plots.jl"))
```

# Miniature microscope

The UCLA Miniscope is a lightweight microscope that utilizes 2-photon fluorescence imaging to record neural activity in awake, freely moving mice [Madruga:2024](@cite). This beginner tutorial aims to reproduce the optical path of the imaging system from the data provided in the [UCLA 2P Miniscope repository](https://github.com/golshanilab/UCLA_2P_Miniscope). You will learn how to:

1. Define optical components (e.g. a [`SphericalDoubletLens`](@ref))
2. Position the mentioned components using the kinematic API
3. Define an optical [`System`](@ref)
4. Add a source of [`Ray`](@ref)s 
5. Trace a [`Beam`](@ref) through the optical system
6. Visualize the results 

All specifications for the optical system and CAD files are taken from this [repository](https://github.com/golshanilab/UCLA_2P_Miniscope) under GNU GPL-3.0.

![Test](test.png)

## Thingy

This tutorial will focus on modeling the return path of the fluorescence light. Based on the design files provided in the resources above, we will start by defining the first lens along the optical path. This plano-convex lens can be generated via three approaches:

1. Design a [`BeamletOptics.AbstractShape`](@ref) and pass it into the [`Lens`](@ref) type
2. Use the [`SphericalLens`](@ref) convenience constructor
3. Utilize the quasi surface-based [`Lens`](@ref) constructor

Due to its versatility, we will use the latter method in order to reimplement the Zemax design file. It is based on the definition of surface types which will be translated into the correct `AbstractShape` automatically. Refer to the [Surface based lens construction](@ref) section for more information. The following code defines the front planar surface and the back spherical convex surface.

**FIXME: explain ref. index**

```julia
# Define the ref. index
NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])

# Objective lens 1
surf_1 = CircularFlatSurface(2*1.144mm)
surf_2 = SphericalSurface(-1.448mm, 2*1.144mm)
obj_lens_1 = Lens(surf_1, surf_2, 1.3mm, NBK7)
```

Note that `obj_lens_1` is a single `Lens` entity once spawned and can be manipulated in 3D-space as described in the section: [Moving optical elements](@ref).

## Building the objective group

To build the second and third lens elements of the objective group we will proceed as above. However, these lenses are spherical doublets formed by two lenses bonded with an optical adhesive. They are also more complex in shape, featuring a mechanical outer diameter. For the second lens the generating code is provided below. Two [`Lens`](@ref)es can be combined into a [`DoubletLens`](@ref). Correct "assembly" of the lens parts is the responsibility of the user.

```julia
# Objective lens 2
surf_1 = SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm)
surf_2 = SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm)
surf_3 = SphericalSurface(-5.020mm, 2*2.380mm)
dl11 = Lens(surf_1, surf_2, 0.5mm, NBK7)
dl12 = Lens(surf_2, surf_3, 2.5mm, NBK7)
translate3d!(dl12, [0, BMO.thickness(dl11), 0])
obj_lens_2 = DoubletLens(dl11, dl12)
```

The lens parts `dl11` and `dl12` are joined together by moving `dl12` along the y-axis by a distance equal to the on-axis thickness of `dl11`.

!!! info "Axis conventions"
    In general, most components provided by **BMO** are aligned with the global y-axis and spawned at the origin unless specified otherwise.

Similarily, the third objective lens can be reproduced from the available Zemax data as follows:

```julia
# Tube lens
surf_1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
surf_2 = SphericalSurface(-3.642mm, 2*3mm)
surf_3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)
dl21 = Lens(surf_1, surf_2, 3.4mm, NBK7)
dl22 = Lens(surf_2, surf_3, 1.0mm, NBK7)
translate3d!(dl22, [0, BMO.thickness(dl21), 0])
tube_lens = DoubletLens(dl21, dl22)
```

Assembling the objective group requires that the individual lens elements be translated into position along the optical axis. This is achieved in the code below by taking the current [`BeamletOptics.position`](@ref) of the elements and adding the on-axis `thickness` in addition to the element distancing taken from the specification.

```julia
translate_to3d!(obj_lens_2, [0, BMO.position(obj_lens_1)[2] + BMO.thickness(obj_lens_1) + 3.344mm, 0])
translate_to3d!(tube_lens, [0, BMO.position(obj_lens_2)[2] + BMO.thickness(obj_lens_2) + 2mm, 0])
objective_group = ObjectGroup([obj_lens_1, obj_lens_2, tube_lens])
```

The `objective_group` locks all elements in place with respect to their relative positions and allows for combined translations and rotations of the elements.
