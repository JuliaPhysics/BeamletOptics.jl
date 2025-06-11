using GLMakie, BeamletOptics

GLMakie.activate!(; ssao=true)

const BMO = BeamletOptics

include(joinpath(@__DIR__, "..", "render_utils.jl"))

const cm = 1e-2
const mm = 1e-3

λ_green = 546.1e-9
λ_red   = 656.3e-9
lambdas = [λ_green, λ_red]
NBK7    = DiscreteRefractiveIndex(lambdas, [1.51872, 1.51432]) # https://www.schott.com/shop/advanced-optics/en/Optical-Glass/SCHOTT-N-BK7/c/glass-SCHOTT%20N-BK7%C2%AE
NSK5    = DiscreteRefractiveIndex(lambdas, [1.59142, 1.58619]) # https://www.schott.com/shop/advanced-optics/en/Optical-Glass/N-SK5/c/glass-N-SK5
NSF4    = DiscreteRefractiveIndex(lambdas, [1.76164, 1.74719]) # https://www.schott.com/shop/advanced-optics/en/Optical-Glass/N-SF4/c/glass-N-SF4
NLAK22  = DiscreteRefractiveIndex(lambdas, [1.65391, 1.64760]) # https://www.schott.com/shop/advanced-optics/en/Optical-Glass/N-LAK22/c/glass-N-LAK22
NLASF44 = DiscreteRefractiveIndex(lambdas, [1.80832, 1.79901]) # https://www.schott.com/shop/advanced-optics/en/Optical-Glass/N-LASF44/c/glass-N-LASF44

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
    obj_lens_1 = Lens(surf_1, surf_2, 1.3mm, NSK5)
    return obj_lens_1
end

function objective_lens_2()
    ## Objective lens 2
    surf_1 = SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm)
    surf_2 = SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm)
    surf_3 = SphericalSurface(-5.020mm, 2*2.380mm)
    dl11 = Lens(surf_1, surf_2, 0.5mm, NSF4)
    dl12 = Lens(surf_2, surf_3, 2.5mm, NLAK22)
    translate3d!(dl12, [0, BMO.thickness(dl11), 0])
    obj_lens_2 = DoubletLens(dl11, dl12)
    return obj_lens_2
end

function tube_doublet()
    ## Tube doublet
    surf_1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
    surf_2 = SphericalSurface(-3.642mm, 2*3mm)
    surf_3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)
    dl21 = Lens(surf_1, surf_2, 3.4mm, NLAK22)
    dl22 = Lens(surf_2, surf_3, 1.0mm, NSF4)
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
    _objective_group = ObjectGroup([obj_lens_1, obj_lens_2, tube_lens])
    xrotate3d!(_objective_group, deg2rad(90))
    return _objective_group
end

function dichroic_splitter()
    ## Dichroic "filter" - glass plate
    shape = BMO.CuboidMesh(8mm, 1mm, 8.5mm)
    translate3d!(shape, [-4mm, 0.0mm, -4.25mm])
    BMO.set_new_origin3d!(shape)
    splitter = Prism(shape, NBK7)
    translate3d!(splitter, [0, 0, 18.677mm])
    xrotate3d!(splitter, deg2rad(45))
    return splitter
end

function emission_filter()
    ef = Prism(BMO.PlanoSurfaceSDF(1mm, 4mm), NBK7)
    return ef
end

function collection_group()
    ## Collection group
    ef_1 = emission_filter()
    ef_2 = emission_filter()

    collect_lens = Lens(
        SphericalSurface(6.580mm, 4.5mm),
        SphericalSurface(-6.580mm, 4.5mm),
        2.6mm,
        NLASF44
    )

    translate3d!(collect_lens, [0, BMO.position(ef_1)[2] + BMO.thickness(ef_1) + 0.1mm, 0])
    translate3d!(ef_2, [0, BMO.position(collect_lens)[2] + BMO.thickness(collect_lens) + 0.25mm, 0])

    collect_group = ObjectGroup([ef_1, collect_lens, ef_2])

    xrotate3d!(collect_group, deg2rad(90))
    translate3d!(collect_group, [0, 0.332mm, 21.937mm])

    return collect_group
end

optical_system_group() = ObjectGroup([objective_group(), dichroic_splitter(), collection_group()])

##
miniscope = miniscope_body()
obj_group = objective_group()
dichro = dichroic_splitter()
col_group = collection_group()

c_view = [
 -0.481504   0.876444  -4.16334e-17   0.00569819
 -0.194339  -0.106766   0.975107     -0.0188647
  0.854627   0.469517   0.221736     -0.0505548
  0.0        0.0        0.0           1.0
]

fig = Figure(size=(600, 600))
display(fig)
brightness = 0.8
ax = LScene(
    fig[1,1],
    scenekw = (
        lights = [
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(-1, 0, 1))
            DirectionalLight(RGBAf(brightness, brightness, brightness), Vec3f(0, -1, -1))
        ], # hella weird syntax but whatever
    )
)
hide_axis(ax)
render!(ax, miniscope, transparency=true, alpha=0.05)
render!(ax, obj_group, transparency=false, color=lens_color())
render!(ax, dichro, transparency=false, color=lens_color())
render!(ax, col_group, transparency=false, color=lens_color())

set_view(ax, c_view)
save("ucla_intro_fig.png", fig; px_per_unit=8, update = false)

## objective lenses
obj_lens_1 = objective_lens_1()
obj_lens_2 = objective_lens_2()
tube_lens = tube_doublet()

## lens 1 render
fig = Figure(size=(600, 200)) 
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, obj_lens_1, transparency=false, color=lens_color())

c_view = [
 -0.704021   0.710179  2.22045e-16  -0.000263177
 -0.412543  -0.408966  0.813975      0.000154901
  0.578068   0.573056  0.5809       -0.00358713
  0.0        0.0       0.0           1.0
]

set_view(ax, c_view)
save("objective_lens_1.png", fig; px_per_unit=8, update = false)

## lens 2 render
fig = Figure(size=(600, 200)) 
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, obj_lens_2, transparency=true, color=lens_color(0.25))

arrow!(ax, [0,0,0], [1,0,0]; color=:red)
arrow!(ax, [0,0,0], [0,1,0]; color=:green)
arrow!(ax, [0,0,0], [0,0,1]; color=:blue)

text!(Point3(5mm, 0, 1mm), text="x")
text!(Point3(0, 5mm, 1mm), text="y")
text!(Point3(0, 0mm, 6.5mm), text="z")

c_view = [
 -0.629489    0.777009  5.55112e-17  -0.000467742
 -0.0282456  -0.022883  0.999339     -0.00229113
  0.776496    0.629073  0.0363517    -0.0134646
  0.0         0.0       0.0           1.0
]

set_view(ax, c_view)
save("objective_lens_2.png", fig; px_per_unit=8, update = false)

## full group render
translate_to3d!(obj_lens_2, [0, BMO.position(obj_lens_1)[2] + BMO.thickness(obj_lens_1) + 3.344mm, 0])
translate_to3d!(tube_lens, [0, BMO.position(obj_lens_2)[2] + BMO.thickness(obj_lens_2) + 2mm, 0])

fig = Figure(size=(600, 300))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, obj_lens_1, transparency=true, color=lens_color(0.5))
render!(ax, obj_lens_2, transparency=true, color=lens_color(0.5))
render!(ax, tube_lens, transparency=true, color=lens_color(0.5))

tp1 = BMO.position(obj_lens_1) + [0, -.5mm, 2mm]
tp2 = BMO.position(obj_lens_2) + [0, .2mm, 3mm]
tp3 = BMO.position(tube_lens) + [0, 1.3mm, 4mm]

text!(tp1, text="obj_lens_1")
text!(tp2, text="obj_lens_2")
text!(tp3, text="tube_lens")

lines!(ax, [0,0], [-1,1], [0,0], color=:black, transparency=true, linestyle=:dashdot)

set_orthographic(ax)

c_view = [
 -0.010335   0.999947      1.02999e-17  -0.00703319
  0.0811125  0.000838342   0.996705     -0.000536583
  0.996651   0.0103009    -0.0811168    -0.460085
  0.0        0.0           0.0           1.0
]

set_view(ax, c_view)
save("full_objective_lens.png", fig; px_per_unit=8, update = false)

## full optical system render
fig = Figure(size=(600,400))
display(fig)
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, obj_group, transparency=false, color=lens_color())
render!(ax, dichro, transparency=false, color=lens_color())
render!(ax, col_group, transparency=false, color=lens_color())

c_view = [
 -0.710179   0.704021  -1.94289e-16  -0.00085023
 -0.261364  -0.26365    0.928535     -0.0135835
  0.653709   0.659426   0.371244     -0.0419803
  0.0        0.0        0.0           1.0
]

set_view(ax, c_view)
save("full_optical_system.png", fig; px_per_unit=8, update = false)

##
osg = optical_system_group()
system = System([osg])
xrotate3d!(osg, deg2rad(-90))

ps_green = PointSource([0, -0.77mm, .25mm], [0, 1, 0], deg2rad(15), λ_green, num_rays=100, num_rings=5)
ps_red = PointSource([0, -0.77mm, -.25mm], [0, 1, 0], deg2rad(15), λ_red, num_rays=100, num_rings=5)

solve_system!(system, ps_green)
solve_system!(system, ps_red)

##
fig = Figure(size=(600,200))
ax = LScene(fig[1,1])
hide_axis(ax)

render!(ax, system; color=lens_color())
render!(ax, ps_green; color=:green2, render_every=1, flen=3mm, show_pos=false)
render!(ax, ps_red; color=:red, render_every=1, flen=3mm, show_pos=false)

display(fig)

set_orthographic(ax)

c_view = [
 0  1   0   -0.015
 0  0   1   0
 1  0   0   -0.55
 0  0   0   1
]

set_view(ax, c_view)

save("miniscope_trace.png", fig; px_per_unit=8, update = false)
