# Aspherical lenses

The package has a basic support for ISO10110 even aspheres. It is planned to extend this support in the future to include extended aspheres and maybe Q-aspheres.
As of now, there is only a convenience constructor for plano-aspheric lenses. However, you can manually construct your own complex shapes and geometries using 
[`SCDI.ConvexAsphericalSurfaceSDF`](@ref) and [`SCDI.ConcaveAsphericalSurfaceSDF`](@ref). Multiple SDFs which represent a closed shape can be added by just adding the SDFs, see [`SCDI.UnionSDF`](@ref).

The following example shows the most simple usage of the plano-aspheric asphere constructor based on the [Thorlabs AL50100J](https://www.thorlabs.com/thorproduct.cfm?partnumber=AL50100J) aspheric lens:

```@example aspheric_lens
using CairoMakie, SCDI

# radius
R = 50.3583e-3
# conic constant
k = -0.789119
# even aspheric coefficients up to 8th order
A = [0, 2.10405e-7*(1e3)^3, 1.76468e-11*(1e3)^5, 1.02641e-15*(1e3)^7]
# center thickness
ct = 10.2e-3
# diameter
d = 50e-3
# refractive index of BK-7 @ 1310 nm (design wavelength)
n = 1.5036

lens = SCDI.PlanoConvexAsphericalLens(R, k, A, d, ct, n)

system = SCDI.System(lens)

fig = Figure(size=(640,480))

ax = fig[1,1] = Axis3(fig, aspect=:data, azimuth=0., elevation=1e-3)
for z in -0.02:0.001:0.02
    pos = [0.0, -0.05, z]
    dir = [0.0, 1.0, 0]
    ray = SCDI.Ray(pos, dir)
    beam = SCDI.Beam(ray)
    SCDI.solve_system!(system, beam, r_max=40)

    SCDI.render_beam!(ax, beam, flen=0.1)
end
SCDI.render_object!(ax, lens)

save("aspherical_lens_showcase.png", fig, px_per_unit=4); nothing # hide

```

![Aspherical lens showcase](aspherical_lens_showcase.png)