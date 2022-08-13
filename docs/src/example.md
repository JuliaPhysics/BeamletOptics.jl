```@example plot
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
using WGLMakie # hide
WGLMakie.activate!() # hide
```

```@example plot

using SCDI

s1 = SCDI.BallLens{Float64}(SCDI.Sphere{Float64}([2., 0, 0], 1), x -> 1.5)
system = SCDI.System([s1])

f = Figure()
ax = f[1, 1] = Axis3(f)
SCDI.render_system!(ax, system)
for z in collect(-0.75:0.1:0.75)
    ray = SCDI.Ray{Float64}([0.0, 0, z], [1.0, 0, 0.0])
    beam = SCDI.Beam([ray], 1e3)
    SCDI.solve_system!(system, beam)
    SCDI.render_beam!(ax, beam, color=:blue, flen=2)
end

f

```