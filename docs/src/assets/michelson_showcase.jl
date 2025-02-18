using GLMakie
using SCDI

cm = 1e-2

splitter_origin = [18.81cm, 23.5cm,0]

asset_dir = joinpath(@__DIR__, "cmi_assets")

NBK7 = SCDI.DiscreteRefractiveIndex([632.8e-9], [1.51509])

## Laser
laser_assembly = SCDI.MeshDummy(joinpath(asset_dir, "Laser Assembly.stl"))

# Mirror
mirror_holder = SCDI.MeshDummy(joinpath(asset_dir, "Mirror Assembly.stl"))
rpm = SCDI.RightAnglePrismMirror(25e-3, 25e-3)
SCDI.zrotate3d!(rpm, deg2rad(45))
SCDI.translate3d!(rpm, [0,33.5cm,0])

mirror_assembly = SCDI.ObjectGroup([rpm, mirror_holder])

# Beamsplitter
splitter_holder = SCDI.MeshDummy(joinpath(asset_dir, "Splitter Assembly.stl"))
cbs = SCDI.CubeBeamsplitter(SCDI.inch, NBK7)
SCDI.zrotate3d!(cbs, deg2rad(-90))

splitter_assembly = SCDI.ObjectGroup([cbs, splitter_holder])

# Arms
arm_holder_1 = SCDI.MeshDummy(joinpath(asset_dir, "Arm Assembly 1.stl"))
arm_holder_2 = SCDI.MeshDummy(joinpath(asset_dir, "Arm Assembly 2.stl"))
m1 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
SCDI.zrotate3d!(m1, deg2rad(-90))
SCDI.translate3d!(m1, [22cm,0,0])
m2 = SCDI.RoundPlanoMirror(SCDI.inch, 5e-3)
SCDI.zrotate3d!(m2, deg2rad(-90))
SCDI.translate3d!(m2, [12cm,0,0])

arm_1 = SCDI.ObjectGroup([m1, arm_holder_1])
arm_2 = SCDI.ObjectGroup([m2, arm_holder_2])

# PD
pd_holder = SCDI.MeshDummy(joinpath(asset_dir, "PD Assembly.stl"))
pd = SCDI.Photodetector(8e-3, 200)
SCDI.translate3d!(pd, [0, -12cm, 0])

pd_assembly = SCDI.ObjectGroup([pd, pd_holder])

system = SCDI.System(
    [
        laser_assembly,
        mirror_assembly,
        splitter_assembly,
        arm_1,
        arm_2,
        pd_assembly
    ]
)

##
SCDI.translate_to3d!(mirror_assembly, [0,-10cm,0])
SCDI.translate_to3d!(splitter_assembly, [18.81cm,23.5cm,0])

SCDI.translate_to3d!(arm_1, splitter_origin)
SCDI.translate_to3d!(arm_2, splitter_origin)
SCDI.translate3d!(arm_1, [3.81cm/2, 0, 0])
SCDI.translate3d!(arm_2, [0, 3.81cm/2, 0])
SCDI.zrotate3d!(arm_2, deg2rad(90))

SCDI.translate_to3d!(pd_assembly, splitter_origin)
SCDI.translate3d!(pd_assembly, [0, -3.81cm/2, 0])

##
fig = Figure()
ax = LScene(fig[1,1])

SCDI.render_object!(ax, laser_assembly)
SCDI.render_object!(ax, mirror_assembly)
SCDI.render_object!(ax, splitter_assembly)
SCDI.render_object!(ax, arm_1)
SCDI.render_object!(ax, arm_2)
SCDI.render_object!(ax, pd_assembly)

beam = SCDI.GaussianBeamlet([0.,0,0], [0.,1,0], 632.8e-9, 5e-4, M2=2)
SCDI.solve_system!(system, beam)

SCDI.render_beam!(ax, beam, color=:red, flen=0.1)

fig

##
i_max = 100
i = 1
p = zeros(i_max)
SCDI.reset_detector!(pd)
while i <= i_max
    SCDI.translate3d!(m1, [5e-9, 0, 0])
    SCDI.solve_system!(system, beam)
    p[i] = SCDI.optical_power(pd)
    SCDI.reset_detector!(pd)
    @info i
    i+=1
end