using GLMakie, BeamletOptics

const BMO = BeamletOptics

GLMakie.activate!(; ssao=true)

Base.include(@__MODULE__, joinpath("..", "render_utils.jl"))

mm = 1e-3
n = 1.5
cbs = RectangularPlateBeamsplitter(36mm, 25mm, 5mm, _->n)
cbs_mount = MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cbs_assembly = ObjectGroup([cbs, cbs_mount])

cmp = RectangularCompensatorPlate(36mm, 25mm, 5mm, _->n)
cmp_mount = MeshDummy(joinpath(file_dir, "PBS Mount.stl"))

cmp_assembly = ObjectGroup([cmp, cmp_mount])

zrotate3d!(cbs_assembly, deg2rad(45))
zrotate3d!(cmp_assembly, deg2rad(-45))

translate3d!(cmp_assembly, [0, 100mm, 0])

system = System([cbs_assembly, cmp_assembly])

beam = GaussianBeamlet([0,-100mm,0], [0, 11, 0], 1e-6, 5e-4)

solve_system!(system, beam)

##
function test(fig, ax)
    splitter_substrate = ax.scene.plots[2]
    splitter_coating = ax.scene.plots[3]
    compensator = ax.scene.plots[5]
    
    splitter_substrate.color[] = :white
    splitter_coating.color[] = :magenta
    compensator.color[] = :white
end

const pbs_view = [
    0.0940287   0.995569   -3.19189e-16  -0.048146
    -0.498165    0.0470503   0.865805      0.0176688
     0.861969   -0.0814105   0.500382     -0.155093
     0.0         0.0         0.0           1.0
]