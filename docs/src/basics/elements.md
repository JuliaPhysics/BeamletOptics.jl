# Optical elements

Optical elements serve as the building blocks for optical systems in the context of this package, representing components such as mirrors, lenses, filters, and more. Unlike the surface/interface based representation of optical elements, they are treated as volumetric bodies in this simulation framework. Optical interactions between rays/beams and elements are defined based on the type of the element and the type of the incident beam/ray. Note that optical elements will also be referred to as [`SCDI.AbstractObject`](@ref)s moving forward.

To ensure compatibility with the [API design](@ref), custom optical elements must adhere to the [`SCDI.AbstractObject`](@ref) interface.

!!! info
    Equations to calculate optical effects often rely on the normal vector at the ray intersection location
    to work correctly. It is important to ensure that this condition is fulfilled when spurious effects occur.

## Types of elements

Some optical elements are provided with this package as is, these include:

- Reflective optical elements
    - [`SCDI.Mirror`](@ref)
- Refractive optical elements
    - [`SCDI.SphericalLens`](@ref)
    - [`SCDI.PlanoConvexAsphericalLens`](@ref)
    - [`SCDI.PlanoConcaveAsphericalLens`](@ref)
    - [`SCDI.Prism`](@ref)
- [`SCDI.Photodetector`](@ref)
- [`SCDI.ThinBeamSplitter`](@ref)

In order to represent geometries of optical elements exactly, this package uses Signed Distance Functions (SDFs) wherever possible and/or feasible. For an introduction into SDFs the [website of Inigo Quilez](https://iquilezles.org/articles/distfunctions/) is referred to. However, other geometry representations are also supported in principle, e.g. [`SCDI.Mesh`](@ref). For more information, refer to the respective documentation. Below, several spherical lenses are showcased:

- [`SCDI.GaussianBeamlet`](@ref) parameters
    - ``w_0 = 5~\text{mm}``
    - ``\lambda=532~\text{nm}``
- Lenses (in order of appearance)
    - [LD1464](https://www.thorlabs.com/thorproduct.cfm?partnumber=LD1464)
    - [LB1811](https://www.thorlabs.com/thorproduct.cfm?partnumber=LB1811)
    - [LC1715](https://www.thorlabs.com/thorproduct.cfm?partnumber=LC1715)
    - [LE1234](https://www.thorlabs.com/thorproduct.cfm?partnumber=LE1234)
    - [LA1805](https://www.thorlabs.com/thorproduct.cfm?partnumber=LA1805)

```@eval
using CairoMakie, SCDI

function NBK7(λ)
    if λ ≈ 532e-9
        return 1.5195
    end
    if λ ≈ 1064e-9
        return 1.5066
    end
    error("Ref. index for λ=$λ not available.")
end

# lens diameter 
d = SCDI.inch

# lens types
r1 = 34.9e-3
r2 = -34.9e-3
l = 6.8e-3
LB1811 = SCDI.SphericalLens(r1, r2, l, d, NBK7)

r1 = Inf
r2 = -15.5e-3
l = 8.6e-3
LA1805 = SCDI.SphericalLens(r1, r2, l, d, NBK7)

r1 = -52e-3
r2 = 52e-3
l = 3e-3
LD1464 = SCDI.SphericalLens(r1, r2, l, d, NBK7)

r1 = Inf
r2 = 25.7e-3
l = 3.5e-3
LC1715 = SCDI.SphericalLens(r1, r2, l, d, NBK7)

r1 = -82.2e-3
r2 = -32.1e-3
l = 3.6e-3
LE1234 = SCDI.SphericalLens(r1, r2, l, d, NBK7)

SCDI.translate3d!(LD1464, [0, 0*d, 0])
SCDI.translate3d!(LB1811, [0, 1*d, 0])
SCDI.translate3d!(LC1715, [0, 2*d, 0])
SCDI.translate3d!(LE1234, [0, 3*d, 0])
SCDI.translate3d!(LA1805, [0, 4*d, 0])

system = SCDI.StaticSystem([
    LB1811,
    LA1805,
    LD1464,
    LC1715,
    LE1234
])

beam = SCDI.GaussianBeamlet(SCDI.Ray([0, -0.05, 0], [0, 1, 0]), 532e-9, 5e-3)
SCDI.solve_system!(system, beam)

fig = Figure(size=(600,240))
aspect = (1,4,1)
limits = (-0.025, 0.025, -0.05, 0.15, -0.025, 0.025)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0., elevation=1e-3)


hidexdecorations!(ax)
hidezdecorations!(ax)

SCDI.render_beam!(ax, beam, color=:green2)
SCDI.render_system!(ax, system)

save("spherical_lens_showcase.png", fig, px_per_unit=4)

nothing
```

The spherical lenses are shown below. To recreate this figure, refer to the [Spherical lenses](@ref) example.

![Spherical lens showcase](spherical_lens_showcase.png)

!!! info "Custom optical elements"
    In order to implement custom geometries and optical elements, refer to the [API design](@ref) section.

## Moving optical elements

Optical elements can move around freely in three-dimensional space, which enables the modeling of kinematics within optical setups. When objects are manipulated, they are translated and rotated around their self-defined center of gravity, which is represented as a ``\mathbb{R}^3``-vector and will be referred to as its [`SCDI.position`](@ref). Additionally, the [`SCDI.orientation`](@ref) of an object, defined as its local fixed coordinate system, is represented by an orthonormal matrix in ``\mathbb{R}^3``. If the object is rotated, this matrix can be used to calculate the inverse transform into global coordinates. 

!!! important
    Elements can be moved freely between each call of [`SCDI.solve_system!`](@ref). However, during tracing it is assumed that all elements remain static.

For elements that implement the [`SCDI.AbstractObject`](@ref) interface, the following movement commands are provided:

- Translation
    - [`SCDI.translate3d!`](@ref)
    - [`SCDI.translate_to3d!`](@ref)
- Rotation
    - [`SCDI.rotate3d!`](@ref)
    - [`SCDI.xrotate3d!`](@ref)
    - [`SCDI.yrotate3d!`](@ref)
    - [`SCDI.zrotate3d!`](@ref)
- Reset commands
    - [`SCDI.reset_translation3d!`](@ref)
    - [`SCDI.reset_rotation3d!`](@ref)

!!! important
    Unless specified otherwise, the translation and rotation commands result in relative motions to the current position and orientation. This must be taken into account when trying to model a specific set of movements.

## Groups of optical elements

For the easier representation of optical elements that move as a group, the [`SCDI.ObjectGroup`](@ref) can be used. Refer to the [Lens groups](@ref) example for more information.

```@docs; canonical=false
SCDI.ObjectGroup
```

