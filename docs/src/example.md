```@example plot
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example plot
using WGLMakie # hide
WGLMakie.activate!() # hide
```

```@example plot

using SCDI, LinearAlgebra

# build a plano-convex lens
normal = [0.5, 0.0, 1.0]
normal /= norm(normal)
scx_lens = SCDI.Lens(
    SCDI.SphericSurface(0.5, 0.25, 0.0),
    SCDI.PlanarSurface(0.25),
    normal,
    [0.0, 0.0, 0.0],
    5e-3,
    x -> 1.5
)

system = SCDI.System([scx_lens])

f = Figure()
ax = f[1, 1] = Axis3(f, azimuth=180, elevation=0)
SCDI.render_system!(ax, system)
for x in collect(-0.2:0.025:0.2)
    ray = SCDI.Ray([x, 0, -0.2], [0.0, 0, 1.0])
    beam = SCDI.Beam([ray], 1e3)
    SCDI.solve_system!(system, beam)
    SCDI.render_beam!(ax, beam, color=:blue, flen=1)
end

xlims!(-.3, .3)
ylims!(-.3, .3)

f

```