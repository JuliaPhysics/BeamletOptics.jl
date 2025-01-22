# Lenses

Lenses are fundamental optical components used to focus or diverge light, making them essential for constructing imaging systems. This package includes a variety of rotationally symmetric lens models to simulate simple imaging setups. All lens models provided as part of this package are based on SDFs. Refer to the [Signed Distance Functions (SDFs)](@ref) section for more information.

## Spherical lenses

Spherical lenses are characterized by surfaces with constant curvature, making them straightforward to model and ideal for basic imaging applications. 

```@docs; canonical=false
SCDI.SphericalLens
```

Below, several spherical lenses are recreated from manufacturer data using the [Spherical lens constructor](@ref):

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

NBK7 = SCDI.DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])

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

The spherical lenses are shown below. To recreate this figure, refer to the [Spherical lens example](@ref).

![Spherical lens showcase](spherical_lens_showcase.png)

### SDF-based spherical lenses

In order to model the lens surfaces shown above, the following SDF-based spherical lens shapes have been implemented:

- [`SCDI.ConvexSphericalSurfaceSDF`](@ref)
- [`SCDI.ConcaveSphericalSurfaceSDF`](@ref)
- [`SCDI.MeniscusLensSDF`](@ref)
- [`SCDI.PlanoSurfaceSDF`](@ref)

They can be combined via the the [`SCDI.UnionSDF`](@ref)-API in order to enable the quasi-surface-based design of spherical lens systems.

!!! hint "Spherical lens example"
    For a complex showcase, refer to the [Double Gauss Lens](@ref) example page.

### Spherical lens constructor

The quasi-surface-based design API mentioned above is based on the [`SCDI.SphericalLensShapeConstructor`](@ref). It builds lenses based on input curvatures, the lens thickness and refractive index as follows:

```@docs; canonical=false
SCDI.SphericalLensShapeConstructor
```

## Doublet lenses

The [`SCDI.DoubletLens`](@ref) is an example for a multi-shape object as mentioned in the [Multi-shape objects](@ref) section.

```@docs; canonical=false
SCDI.DoubletLens
```

The following image shows the [AC254-150-AB](https://www.thorlabs.com/thorproduct.cfm?partnumber=AC254-150-AB) doublet lens for 488 and 707 nm. It has been created using the [`SCDI.SphericalDoubletLens`](@ref) constructor. Changes in the [`SCDI.Ray`](@ref) direction are marked with red dots (disregarding the beam spawn point).


```@eval
using CairoMakie, SCDI

λs = [488e-9, 707e-9, 1064e-9]

NLAK22 = SCDI.DiscreteRefractiveIndex(λs, [1.6591, 1.6456, 1.6374])
NSF10 = SCDI.DiscreteRefractiveIndex(λs, [1.7460, 1.7168, 1.7021])

AC254_150_AB = SCDI.SphericalDoubletLens(87.9e-3, 105.6e-3, 1000, 6e-3, 3e-3, SCDI.inch, NLAK22, NSF10)

system = SCDI.System([AC254_150_AB])

fig = Figure(size=(600,240))
aspect = (1,4,1)
limits = (-0.025, 0.025, -0.025, 0.175, -0.025, 0.025)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0., elevation=1e-3)

hidexdecorations!(ax)
hidezdecorations!(ax)

SCDI.render_system!(ax, system)

zs_1 = LinRange(-0.011, 0.011, 6)
zs_2 = LinRange(-0.01, 0.01, 5)

for (i, z) in enumerate(zs_1)
    beam = SCDI.Beam([0, -0.02 , z], [0,1.,0], 488e-9)
    SCDI.solve_system!(system, beam)
    SCDI.render_beam!(ax, beam, flen=0.15, color=RGBAf(0,0,1,0.7))
end

for (i, z) in enumerate(zs_2)
    beam = SCDI.Beam([0, -0.02 , z], [0,1.,0], 707e-9)
    SCDI.solve_system!(system, beam)
    SCDI.render_beam!(ax, beam, flen=0.15, color=RGBAf(1,0,0,0.5), show_pos=true)
end

save("doublet_showcase.png", fig, px_per_unit=4)

nothing
```

![Doublet lens showcase](doublet_showcase.png)

## Aspherical lenses

Aspherical lenses offer more advanced control over aberrations, enabling higher performance in specialized optical systems. The package offers surface support for rotationally symetrical [aspheric lenses](https://en.wikipedia.org/wiki/Aspheric_lens) that adhere to the DIN ISO 10110 convention with even terms. This includes the following types:

- [`SCDI.PlanoConcaveAsphericalLens`](@ref)
- [`SCDI.PlanoConvexAsphericalLens`](@ref)

!!! hint "Aspherical lens example"
    Refer to the [Aspherical lens example](@ref) for a showcase on how to implement a plano-convex asphere.

```@docs; canonical=false
SCDI.PlanoConvexAsphericalLens
```