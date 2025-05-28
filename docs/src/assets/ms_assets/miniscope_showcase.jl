using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

const cm = 1e-2
const mm = 1e-3

const NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])

##
function miniscope_body()
    miniscope = MeshDummy(joinpath(@__DIR__, "Miniscope.stl"))
    # offset from CAD model
    dx = -3.646mm 
    dy = 41.933mm
    dz = -7.911mm
    # fix mesh offset
    translate_to3d!(miniscope, [dx, dy, dz])
    BMO.set_new_origin3d!(miniscope)
    # fix mesh rotation
    xrotate3d!(miniscope, deg2rad(90))
    zrotate3d!(miniscope, deg2rad(180))
    BMO.set_new_origin3d!(miniscope)
    return miniscope
end

function objective_lens_1()
    ## Objective lens 1
    surf_1 = CircularFlatSurface(2*1.144mm)
    surf_2 = SphericalSurface(-1.448mm, 2*1.144mm)
    obj_lens_1 = Lens(surf_1, surf_2, 1.3mm, NBK7)
    return obj_lens_1
end

function objective_lens_2()
    ## Objective lens 2
    surf_1 = SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm)
    surf_2 = SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm)
    surf_3 = SphericalSurface(-5.020mm, 2*2.380mm)
    dl11 = Lens(surf_1, surf_2, 0.5mm, NBK7)
    dl12 = Lens(surf_2, surf_3, 2.5mm, NBK7)
    translate3d!(dl12, [0, BMO.thickness(dl11), 0])
    obj_lens_2 = DoubletLens(dl11, dl12)
    return obj_lens_2
end

function tube_doublet()
    ## Tube doublet
    surf_1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
    surf_2 = SphericalSurface(-3.642mm, 2*3mm)
    surf_3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)
    dl21 = Lens(surf_1, surf_2, 3.4mm, NBK7)
    dl22 = Lens(surf_2, surf_3, 1.0mm, NBK7)
    translate3d!(dl22, [0, BMO.thickness(dl21), 0])
    tube_lens = DoubletLens(dl21, dl22)
    return tube_lens
end

function objective_group()
    obj_lens_1 = objective_lens_1()
    obj_lens_2 = objective_lens_2()
    tube_lens = tube_doublet()
    translate_to3d!(obj_lens_2, [0, BMO.position(obj_lens_1)[2] + BMO.thickness(obj_lens_1) + 3.344mm, 0])
    translate_to3d!(tube_lens, [0, BMO.position(obj_lens_2)[2] + BMO.thickness(obj_lens_2) + 2mm, 0])
    objective_group = ObjectGroup([obj_lens_1, obj_lens_2, tube_lens])
    xrotate3d!(objective_group, deg2rad(90))
    return objective_group
end

function dichroic_filter()
    ## Dichroic "filter" - glass plate
    shape = BMO.CuboidMesh(8mm, 1mm, 8.5mm)
    translate3d!(shape, [-4mm, 0.0mm, -4.25mm])
    BMO.set_new_origin3d!(shape)
    translate3d!(shape, [0, 0, 18.677mm])
    xrotate3d!(shape, deg2rad(45))
    filter = Lens(shape, NBK7)
    return filter
end

##
miniscope = miniscope_body()
filter = dichroic_filter()

c_view = [
 -0.481504   0.876444  -4.16334e-17   0.00569819
 -0.194339  -0.106766   0.975107     -0.0188647
  0.854627   0.469517   0.221736     -0.0505548
  0.0        0.0        0.0           1.0
]

fig = Figure(size=(600, 600))
display(fig)
brightness = .8
ax = LScene(
    fig[1,1],
    scenekw = (
        lights = [
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(1, 0, 1))
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(0, 1, 1))
        ], # hella weird syntax but whatever
    )
)
hide_axis(ax)
render!(ax, miniscope)
render!(ax, filter; color=RGBAf(0,0,1,0.25))

set_view(ax, c_view)
save("test.png", fig; px_per_unit=8, update = false)

##
lg = objective_group()

xrotate3d!(lg, deg2rad(90))

l1 = lg.objects[1]
l2 = lg.objects[2]
l3 = lg.objects[3]


fig = Figure(size=(600, 300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, l1, color=RGBAf(1, 0, 0, .25))
render!(ax, l2, color=RGBAf(0, 1, 0, .25))
render!(ax, l3, color=RGBAf(0, 0, 1, .25))

set_orthographic(ax)

c_view = [
  0.00382847  -0.999993    1.11998e-16  -0.00708868
  0.118809     0.00045486  0.992917      4.83965e-5
 -0.99291     -0.00380136  0.11881      -0.432527
  0.0          0.0         0.0           1.0
]

set_view(ax, c_view)