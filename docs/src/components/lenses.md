# Lenses

Lenses are fundamental optical components used to focus or diverge light, making them essential for constructing imaging systems. The [`BeamletOptics.AbstractRefractiveOptic`](@ref) type provides a general definition of components that refract light. This package includes a variety of rotationally symmetric lens models to simulate simple imaging setups. All lens models provided as part of this package are based on SDFs. Refer to the [Signed Distance Functions (SDFs)](@ref) section for more information.

## Spherical lenses

Spherical lenses are characterized by surfaces with constant curvature, making them straightforward to model and ideal for basic imaging applications. 

```@docs; canonical=false
SphericalLens
```

Below, several spherical lenses are recreated from manufacturer data using the [Spherical lens constructor](@ref):

- [`GaussianBeamlet`](@ref) parameters
    - ``w_0 = 5~\text{mm}``
    - ``\lambda=532~\text{nm}``
- Lenses (in order of appearance)
    - [LD1464](https://www.thorlabs.com/thorproduct.cfm?partnumber=LD1464)
    - [LB1811](https://www.thorlabs.com/thorproduct.cfm?partnumber=LB1811)
    - [LC1715](https://www.thorlabs.com/thorproduct.cfm?partnumber=LC1715)
    - [LE1234](https://www.thorlabs.com/thorproduct.cfm?partnumber=LE1234)
    - [LA1805](https://www.thorlabs.com/thorproduct.cfm?partnumber=LA1805)

The spherical lenses are shown below. To recreate this figure, refer to the [Spherical lens example](@ref).

```@eval
file_dir = joinpath(@__DIR__, "..", "assets")

Base.include(@__MODULE__, joinpath(file_dir, "spherical_lens_showcase.jl"))

save("spherical_lens_showcase.png", fig, px_per_unit=4); nothing
```

![Spherical lens showcase](spherical_lens_showcase.png)

### SDF-based spherical lenses

In order to model the lens surfaces shown above, the following SDF-based spherical lens shapes have been implemented:

- [`BeamletOptics.ConvexSphericalSurfaceSDF`](@ref)
- [`BeamletOptics.ConcaveSphericalSurfaceSDF`](@ref)
- [`BeamletOptics.MeniscusLensSDF`](@ref)
- [`BeamletOptics.PlanoSurfaceSDF`](@ref)

In general, [`BeamletOptics.AbstractLensSDF`](@ref)s can be combined via the the [`BeamletOptics.UnionSDF`](@ref)-API in order to enable the quasi-surface-based design of spherical lens systems.

!!! hint "Spherical lens example"
    For a complex showcase, refer to the [Double Gauss lens](@ref) example page.

### Spherical lens constructor

The quasi-surface-based design API mentioned above is based on the [`BeamletOptics.SphericalLensShapeConstructor`](@ref). It builds lenses based on input curvatures, the lens thickness and refractive index as follows:

```@docs; canonical=false
BeamletOptics.SphericalLensShapeConstructor
```

## Doublet lenses

The [`DoubletLens`](@ref) is an example for a multi-shape object as mentioned in the [Multi-shape objects](@ref) section.

```@docs; canonical=false
SphericalDoubletLens(::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any, ::Any)
```

The following image shows the [AC254-150-AB](https://www.thorlabs.com/thorproduct.cfm?partnumber=AC254-150-AB) doublet lens for 488 and 707 nm. It has been created using the [`SphericalDoubletLens`](@ref) constructor shown above. Changes in the [`Ray`](@ref) direction are marked with red dots (disregarding the beam spawn point).


```@eval
using CairoMakie, BeamletOptics

λs = [488e-9, 707e-9, 1064e-9]

NLAK22 = BeamletOptics.DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
NSF10 = BeamletOptics.DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

AC254_150_AB = SphericalDoubletLens(87.9e-3, 105.6e-3, 1000, 6e-3, 3e-3, BeamletOptics.inch, NLAK22, NSF10)

system = System([AC254_150_AB])

fig = Figure(size=(600,240))
aspect = (1,4,1)
limits = (-0.025, 0.025, -0.025, 0.175, -0.025, 0.025)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0., elevation=1e-3)

hidedecorations!(ax)
hidespines!(ax)

render_system!(ax, system)

zs_1 = LinRange(-0.011, 0.011, 6)
zs_2 = LinRange(-0.01, 0.01, 5)

for (i, z) in enumerate(zs_1)
    beam = Beam([0, -0.02 , z], [0,1.,0], 488e-9)
    solve_system!(system, beam)
    render_beam!(ax, beam, flen=0.15, color=RGBAf(0,0,1,0.7))
end

for (i, z) in enumerate(zs_2)
    beam = Beam([0, -0.02 , z], [0,1.,0], 707e-9)
    solve_system!(system, beam)
    render_beam!(ax, beam, flen=0.15, color=RGBAf(1,0,0,0.5), show_pos=true)
end

save("doublet_showcase.png", fig, px_per_unit=4)

nothing
```

![Doublet lens showcase](doublet_showcase.png)

## Aspherical lenses

Aspherical lenses offer more advanced control over aberrations, enabling higher performance in specialized optical systems. The package offers surface support for rotationally symetrical [aspheric lenses](https://en.wikipedia.org/wiki/Aspheric_lens) that adhere to the DIN ISO 10110 convention with even terms.

To construct a lens with any possible combination of convex/concave, spherical/aspherical surfaces you can use the `Lens` constructor. A complex example of such a lens might look like the following example. This lens has the following peculiarities:
- The front surface is an aspherical convex surface with a clear diameter smaller than the full mechanical diameter
- The back surface is an aspherical concave surface which first curves outwards before change slope and curving invards, giving a more "convex" like character while still beeing a concave lens by definition. Also this surface extends towards the full outer diameter.

```@example
using CairoMakie, BeamletOptics

L3 = Lens(
        EvenAsphericSurface(
            3.618e-3, # r
            3.04e-3, # d
            -44.874, # conic
            [0,-0.14756*(1e3)^3, 0.035194*(1e3)^5, -0.0032262*(1e3)^7,
            0.0018592*(1e3)^9, 0.00036658*(1e3)^11, -0.00016039*(1e3)^13,
            -3.1846e-5*(1e3)^15] # coeffs
        ),
        EvenAsphericSurface(
            2.161e-3, # r
            3.7e-3, # d
            -10.719, # conic
            [0,-0.096568*(1e3)^3, 0.026771*(1e3)^5, -0.011261*(1e3)^7,
            0.0019879*(1e3)^9, 0.00015579*(1e3)^11, -0.00012433*(1e3)^13,
            1.5264e-5*(1e3)^15] # coeffs
        ),
        0.7e-3, # center_thickness
        n -> 1.580200
    )

system = System([L3])

fig = Figure(size=(600,240))
ax = Axis3(fig[1,1], aspect=:data, azimuth=0., elevation=1e-3)

hidedecorations!(ax)
hidespines!(ax)

render_system!(ax, system)

fig
```

!!! hint "Aspherical lens example"
    Refer to the [Simple aspherical lens example](@ref) for a showcase on how to implement a plano-convex asphere.

```@docs; canonical=false
BeamletOptics.Lens
```