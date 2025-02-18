using CairoMakie, SCDI

cbs = SCDI.CubeBeamsplitter(10e-3, n->1.5)
SCDI.zrotate3d!(cbs, deg2rad(90 + 15))

fig = Figure(size=(600,350))
aspect = (1,2,2)
limits = (-0.01+2.5e-3, 0.01+2.5e-3, -0.02, 0.02, -0.02, 0.02)
ax = Axis3(fig[1,1], aspect=aspect, limits=limits, azimuth=0., elevation=pi/2)

hidedecorations!(ax)
hidespines!(ax)

SCDI.render_object!(ax, cbs)

beam = SCDI.GaussianBeamlet([0,-20e-3,0], [0,1,0], 400e-9, .2e-3)
system = SCDI.System([cbs])

SCDI.solve_system!(system, beam)

SCDI.render_beam!(ax, beam, color=:blue, flen=20e-3)
