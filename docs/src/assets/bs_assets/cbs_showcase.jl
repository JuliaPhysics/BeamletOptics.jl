using GLMakie, SCDI

file_dir = @__DIR__

mm = 1e-3
n = 1.5

cbs1 = SCDI.CubeBeamsplitter(SCDI.inch, _->n)
cbs1_mount = SCDI.MeshDummy(joinpath(file_dir, "CBS Mount.stl"))

cbs2 = SCDI.CubeBeamsplitter(SCDI.inch, _->n)
cbs2_mount = SCDI.MeshDummy(joinpath(file_dir, "CBS Mount.stl"))
SCDI.zrotate3d!(cbs2_mount, π)

cbs1_assembly = SCDI.ObjectGroup([cbs1, cbs1_mount])
cbs2_assembly = SCDI.ObjectGroup([cbs2, cbs2_mount])

m1 = SCDI.RightAnglePrismMirror(SCDI.inch, SCDI.inch)
m1_mount = SCDI.MeshDummy(joinpath(file_dir, "CBS Mount.stl"))
SCDI.zrotate3d!(m1, π)

m2 = SCDI.RightAnglePrismMirror(SCDI.inch, SCDI.inch)
m2_mount = SCDI.MeshDummy(joinpath(file_dir, "CBS Mount.stl"))
SCDI.zrotate3d!(m2, π)

m1_assembly = SCDI.ObjectGroup([m1, m1_mount])
m2_assembly = SCDI.ObjectGroup([m2, m2_mount])

SCDI.translate_to3d!(cbs2_assembly, [-100mm, 150mm, 0])

SCDI.translate_to3d!(m1_assembly, [-100mm, 0, 0])
SCDI.zrotate3d!(m1_assembly, deg2rad(-45))

SCDI.translate_to3d!(m2_assembly, [0, 150mm, 0])
SCDI.zrotate3d!(m2_assembly, deg2rad(135))

system = SCDI.System([cbs1_assembly, cbs2_assembly, m1_assembly, m2_assembly])

beam = SCDI.GaussianBeamlet([0,-60mm,0], [0, 1, 0], 1e-6, 1e-3)

SCDI.solve_system!(system, beam)

fig = Figure(size=(600, 470))
ax = Axis3(fig[1,1]; aspect=:data,  azimuth=0, elevation=pi/2)

SCDI.render_object!.(ax, [cbs1_mount, cbs2_mount, m1_mount, m2_mount])
SCDI.render_object!.(ax, [cbs1, cbs2, m1, m2])
SCDI.render_beam!(ax, beam, flen=40mm, color=:red)

hidedecorations!(ax)
hidespines!(ax)