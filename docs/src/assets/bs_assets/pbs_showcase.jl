using CairoMakie, SCDI

mm = 1e-3
n = 1.5
cbs = SCDI.RectangularPlateBeamsplitter(36mm, 25mm, 5mm, _->n)
cbs_mount = SCDI.MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cbs_assembly = SCDI.ObjectGroup([cbs, cbs_mount])

cmp = SCDI.RectangularCompensatorPlate(36mm, 25mm, 5mm, _->n)
cmp_mount = SCDI.MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cmp_assembly = SCDI.ObjectGroup([cmp, cmp_mount])

SCDI.zrotate3d!(cbs_assembly, deg2rad(-45))
SCDI.zrotate3d!(cmp_assembly, deg2rad(45))

SCDI.translate3d!(cmp_assembly, [0, 100mm, 0])

system = SCDI.System([cbs_assembly, cmp_assembly])

beam = SCDI.GaussianBeamlet([0,-40mm,0], [0, 11, 0], 1e-6, 5e-4)

SCDI.solve_system!(system, beam)

fig = Figure(size=(600, 320))
ax = Axis3(fig[1,1]; aspect=:data,  azimuth=0, elevation=pi/2)

SCDI.render_object!(ax, cbs_mount)
SCDI.render_object!(ax, cbs)
SCDI.render_object!(ax, cmp_mount)
SCDI.render_object!(ax, cmp)
SCDI.render_beam!(ax, beam, flen=60mm, color=:red)

hidedecorations!(ax)
hidespines!(ax)

fp = ax.scene.plots[14]
ct = ax.scene.plots[15]
cm = ax.scene.plots[17]

fp.color[] = :white
ct.color[] = :magenta
cm.color[] = :white

fig