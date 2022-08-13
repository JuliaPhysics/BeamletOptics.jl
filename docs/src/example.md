```@example plot
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
using WGLMakie # hide
WGLMakie.activate!() # hide
```

```@example plot

using SCDI

struct Plane{T} <: SCDI.AbstractMesh
    mesh::SCDI.Mesh{T}
end

vertices = [
        1 1 0
        1 -1 0
        -1 -1 0
        -1 1 0
    ]
    faces = [
        1 2 3
        3 4 1
    ]
    pos = [0, 0, 0]
    dir = [1 0 0;0 1 0; 0 0 1]
    scale = 1
    plane = Plane{Float64}(SCDI.Mesh{Float64}(
        vertices,
        faces,
        dir,
        pos,
        scale
    ))

system = SCDI.System([plane])

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