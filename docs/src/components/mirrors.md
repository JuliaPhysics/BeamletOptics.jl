# Mirrors

## Plano Mirrors

```@eval
using CairoMakie, SCDI

asset_dir = joinpath(@__DIR__, "../assets")

function spawn_mirror_mount()
    m1 = SCDI.RoundPlanoMirror(1SCDI.inch, 6e-3)
    holder = SCDI.MeshDummy(joinpath(asset_dir, "Mirror_Post.stl"))
    SCDI.translate_to3d!(holder, [0,0,-(5.68e-2)])
    SCDI.set_new_origin3d!(holder)
    mirror_mount = SCDI.ObjectGroup([holder, m1])
    return mirror_mount
end

m1 = spawn_mirror_mount()
m2 = spawn_mirror_mount()
m3 = spawn_mirror_mount()
m4 = spawn_mirror_mount()

SCDI.zrotate3d!(m1, deg2rad(45))
SCDI.zrotate3d!(m2, deg2rad(-135))
SCDI.zrotate3d!(m3, deg2rad(-45))
SCDI.zrotate3d!(m4, deg2rad(135))

dy = 0.3

SCDI.translate_to3d!(m1, [0,    0,      0])
SCDI.translate_to3d!(m2, [0.1,  0,      0])
SCDI.translate_to3d!(m3, [0.1,  dy,     0])
SCDI.translate_to3d!(m4, [0,    dy,     0])

##
system = SCDI.StaticSystem([m1, m2, m3, m4])

fig = Figure(size=(600, 300))
limits = (-0.1, 0.2, -0.1, 0.4, -0.075, 0)
aspect = (3,5,0.75)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0.0, elevation=π/4)

hidedecorations!(ax)

beam = SCDI.GaussianBeamlet(SCDI.Ray([0, -0.1, 0], [0, 1, 0]), 532e-9, 1e-3)

SCDI.solve_system!(system, beam, r_max=100)

SCDI.render_object!(ax, m1)
SCDI.render_object!(ax, m4)
SCDI.render_beam!(ax, beam, flen=0.1, color=RGBf(0,1,0))
SCDI.render_object!(ax, m2)
SCDI.render_object!(ax, m3)

save("plano_mirror_showcase.png", fig, px_per_unit=4); nothing
```

![Plano mirror showcase](plano_mirror_showcase.png)

## Concave Mirrors

Spherical [`SCDI.ConcaveSphericalMirror`](@ref)s can be used to focus light via reflection. The 

```@eval 
using CairoMakie, SCDI

# new figure

distance = 20e-2
factor = 1.2
RoC = distance/2 * factor
m1 = SCDI.ConcaveSphericalMirror(RoC, 5e-3, 2SCDI.inch)
m2 = SCDI.ConcaveSphericalMirror(RoC, 5e-3, 2SCDI.inch)

SCDI.zrotate3d!(m1, deg2rad(180))
SCDI.translate3d!(m2, [0, distance, 0])

system = SCDI.StaticSystem([m1, m2])

fig = Figure(size=(600,240))
dr = 0.03
y1 = -0.02
y2 = 0.22
aspect = (1, (y2-y1)/(2dr), 1)
limits = (-dr, dr, y1, y2, -dr, dr)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0.4, elevation=1e-3)

hidedecorations!(ax)

beam = SCDI.Beam(SCDI.Ray([0, distance/2, 7e-3], [0.17, 1, 0]))

SCDI.solve_system!(system, beam, r_max=100)

SCDI.render_object!(ax, m1)
SCDI.render_beam!(ax, beam, flen=0.1)
SCDI.render_object!(ax, m2)

save("concave_mirror_showcase.png", fig, px_per_unit=4); nothing
```

![Concave mirror multipass showcase](concave_mirror_showcase.png)