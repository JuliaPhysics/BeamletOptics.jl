include(joinpath(@__DIR__, "..", "render_utils.jl"))

function reset_beamlet!(beam::GaussianBeamlet)
    BeamletOptics._drop_beams!(beam)
    r_c = first(BeamletOptics.rays(beam.chief))
    r_w = first(BeamletOptics.rays(beam.waist))
    r_d = first(BeamletOptics.rays(beam.divergence))
    BeamletOptics.intersection!(r_c, nothing)
    BeamletOptics.intersection!(r_w, nothing)
    BeamletOptics.intersection!(r_d, nothing)
    beam.chief.rays = [r_c]
    beam.waist.rays = [r_w]
    beam.divergence.rays = [r_d]
    return nothing
end

const cm = 1e-2
const mm = 1e-3

splitter_origin = [18.81cm, 23.5cm,0]

asset_dir = @__DIR__

NBK7 = DiscreteRefractiveIndex([632.8e-9], [1.51509])

## Laser
laser_assembly = MeshDummy(joinpath(asset_dir, "Laser Assembly.stl"))

# Mirror
mirror_holder = MeshDummy(joinpath(asset_dir, "Mirror Assembly.stl"))
rpm = RightAnglePrismMirror(25e-3, 25e-3)
zrotate3d!(rpm, deg2rad(45))
translate3d!(rpm, [0,33.5cm,0])

mirror_assembly = ObjectGroup([rpm, mirror_holder])

# Beamsplitter
splitter_holder = MeshDummy(joinpath(asset_dir, "Splitter Assembly.stl"))
cbs = CubeBeamsplitter(BeamletOptics.inch, NBK7)
zrotate3d!(cbs, deg2rad(-90))

splitter_assembly = ObjectGroup([cbs, splitter_holder])

# Arms
arm_holder_1 = MeshDummy(joinpath(asset_dir, "Arm Assembly 1.stl"))
arm_holder_2 = MeshDummy(joinpath(asset_dir, "Arm Assembly 2.stl"))
m1 = RoundPlanoMirror(BeamletOptics.inch, 5e-3)
zrotate3d!(m1, deg2rad(-90))
translate3d!(m1, [22cm,0,0])
m2 = RoundPlanoMirror(BeamletOptics.inch, 5e-3)
zrotate3d!(m2, deg2rad(-90))
translate3d!(m2, [12cm,0,0])

arm_1 = ObjectGroup([m1, arm_holder_1])
arm_2 = ObjectGroup([m2, arm_holder_2])

# PD
pd_holder = MeshDummy(joinpath(asset_dir, "PD Assembly.stl"))
pd = Photodetector(8e-3, 200)
translate3d!(pd, [0, -12cm, 0])

pd_assembly = ObjectGroup([pd, pd_holder])

system = System(
    [
        laser_assembly,
        mirror_assembly,
        splitter_assembly,
        arm_1,
        arm_2,
        pd_assembly
    ]
)

## setup system (dummies and optics)
translate_to3d!(mirror_assembly, [0,-10cm,0])
translate_to3d!(splitter_assembly, [18.81cm,23.5cm,0])

translate_to3d!(arm_1, splitter_origin)
translate_to3d!(arm_2, splitter_origin)
translate3d!(arm_1, [3.81cm/2, 0, 0])
translate3d!(arm_2, [0, 3.81cm/2, 0])
zrotate3d!(arm_2, deg2rad(90))

translate_to3d!(pd_assembly, splitter_origin)
translate3d!(pd_assembly, [0, -3.81cm/2, 0])