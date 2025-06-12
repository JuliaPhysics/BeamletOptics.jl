using GLMakie, BeamletOptics

const BMO = BeamletOptics

NBK7 = DiscreteRefractiveIndex([532e-9, 1064e-9], [1.5195, 1.5066])

const mm = 1e-3

hide_axis(ls::LScene) = (ls.show_axis[] = false)

##
miniscope = MeshDummy(joinpath(@__DIR__, "Miniscope.stl"))

dx = -3.646mm 
dy = 41.933mm
dz = -7.911mm

translate_to3d!(miniscope, [dx, dy, dz])

BMO.set_new_origin3d!(miniscope)

xrotate3d!(miniscope, deg2rad(90))
zrotate3d!(miniscope, deg2rad(180))

BMO.set_new_origin3d!(miniscope)

## Objective lens 1
surf_1 = CircularFlatSurface(2*1.144mm)
surf_2 = SphericalSurface(-1.448mm, 2*1.144mm)

obj_lens_1 = Lens(surf_1, surf_2, 1.3mm, NBK7)

## Objective lens 2
surf_1 = SphericalSurface(38.184mm, 2*1.840mm, 2*2.380mm)
surf_2 = SphericalSurface(3.467mm, 2*2.060mm, 2*2.380mm)
surf_3 = SphericalSurface(-5.020mm, 2*2.380mm)

dl11 = Lens(surf_1, surf_2, 0.5mm, NBK7)
dl12 = Lens(surf_2, surf_3, 2.5mm, NBK7)

translate3d!(dl12, [0, thickness(dl11), 0])

obj_lens_2 = DoubletLens(dl11, dl12)

## Tube doublet
mech_diam = true

if mech_diam
    surf_1 = SphericalSurface(7.744mm, 2*2.812mm, 2*3mm)
    surf_2 = SphericalSurface(-3.642mm, 2*3mm)
    surf_3 = SphericalSurface(-14.413mm, 2*2.812mm, 2*3mm)
else
    surf_1 = SphericalSurface(7.744mm, 2*3mm)
    surf_2 = SphericalSurface(-3.642mm, 2*3mm)
    surf_3 = SphericalSurface(-14.413mm, 2*3mm)
end

dl21 = Lens(surf_1, surf_2, 3.4mm, NBK7)
dl22 = Lens(surf_2, surf_3, 1.0mm, NBK7)

translate3d!(dl22, [0, thickness(dl21), 0])

tube_lens = DoubletLens(dl21, dl22)

## Objective group
translate_to3d!(obj_lens_2, [0, BMO.position(obj_lens_1)[2] + thickness(obj_lens_1) + 3.344mm, 0])
translate_to3d!(tube_lens, [0, BMO.position(obj_lens_2)[2] + thickness(obj_lens_2) + 2mm, 0])

objective_group = ObjectGroup([obj_lens_1, obj_lens_2, tube_lens])

xrotate3d!(objective_group, deg2rad(90))

## Dichroic "filter"
shape = BMO.CuboidMesh(8mm, 1mm, 8.5mm)
translate3d!(shape, [-4mm, 0.0mm, -4.25mm])
BMO.set_new_origin3d!(shape)

translate3d!(shape, [0, 0, 18.677mm])
xrotate3d!(shape, deg2rad(45))

filter = Lens(shape, NBK7)

## Collection group
ef_1 = Lens(BMO.PlanoSurfaceSDF(1mm, 4mm), NBK7)
ef_2 = Lens(BMO.PlanoSurfaceSDF(1mm, 4mm), NBK7)

collect_lens = Lens(
    SphericalSurface(6.580mm, 4.5mm),
    SphericalSurface(-6.580mm, 4.5mm),
    2.6mm,
    NBK7
)

translate3d!(collect_lens, [0, BMO.position(ef_1)[2] + thickness(ef_1) + 0.1mm, 0])
translate3d!(ef_2, [0, BMO.position(collect_lens)[2] + thickness(collect_lens) + 0.25mm, 0])

collect_group = ObjectGroup([ef_1, collect_lens, ef_2])

xrotate3d!(collect_group, deg2rad(90))
translate3d!(collect_group, [0, 0.332mm, 21.937mm])

##
fig = Figure()
ax = LScene(fig[1,1])
hide_axis(ax)
render!(ax, miniscope, transparency=true, alpha=0.05)

system = System([objective_group, filter, collect_group])

render!(ax, system)

ys = LinRange(-.3, .3, 30)

##
for y in ys
    beam = Beam([0, 0, -0.77mm], [0, y, 1], 1.064e-6)

    solve_system!(system, beam)
    render!(ax, beam, flen=4mm, color=:blue)
end

for y in ys
    beam = Beam([0, .5mm, -0.77mm], [0, y, 1], 1.064e-6)

    solve_system!(system, beam)
    render!(ax, beam, flen=4mm, show_pos=false, color=:red)
end

for y in ys
    beam = Beam([0, -.5mm, -0.77mm], [0, y, 1], 1.064e-6)

    solve_system!(system, beam)
    render!(ax, beam, flen=4mm, show_pos=false, color=:green)
end

fig
