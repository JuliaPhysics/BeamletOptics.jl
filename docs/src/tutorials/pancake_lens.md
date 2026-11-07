```@meta
EditURL = "../../literate/pancake_lens.jl"
```

# Pancake lens

!!! tip "Before you start"
    Install the packages of this tutorial once: press `]` in the REPL and type
    `add BeamletOptics, GLMakie, AbstractTrees`. The [Laser alignment](@ref) tutorial introduces the
    basics of placing components and tracing beams.

The eyepiece of a compact VR headset is often a *pancake lens*: the light is folded three times
within a lens a few millimeters thick, so the display sits close to the eye. The fold is a
polarization trick that only works because several thin films act on one surface, in a fixed order,
and the light passes them more than once. This tutorial builds such an eyepiece from a single
plano-convex lens and traces it non-sequentially: the solver finds the folded path by itself.

```@raw html
<div class="bmo-card">
<p><Badge type="danger" text="Advanced" /></p>
```

You will learn how to:

1. Coat the faces of a lens with [`surface_model!`](@ref)
2. Stack polarizing films with [`SurfaceStack`](@ref), [`Retarder`](@ref) and [`ReflectivePolarizer`](@ref)
3. Trace beam trees that split at partial mirrors, and sort their branches
4. Find the ghost image of a pancake lens, and how it depends on retarder errors and color

```@raw html
</div>
```

```@raw html
<div class="bmo-card">
<p><a href="../pancake_lens.jl">Download tutorial <code>.jl</code></a> &middot; <a href="../pancake_lens.ipynb">Download Jupyter notebook</a></p>
</div>
```

## How the fold works

The lens carries two coatings. Its convex front face, towards the display, is a 50:50 mirror. Its flat
back face, towards the eye, carries a quarter-wave retarder film under a reflective polarizer, which
transmits light polarized along x and reflects the perpendicular component. The display emits
circularly polarized light.

1. Half of the light passes the 50:50 mirror into the lens.
2. The retarder turns the circular light into linear light along z, which the polarizer **reflects**.
3. Back through the retarder the light is circular again. The 50:50 mirror reflects half of it and
   flips its handedness; seen from inside the glass the mirror is concave and gives the eyepiece most
   of its optical power.
4. The retarder now turns the light into linear light along x, which the polarizer **transmits**.

At best a quarter of the light reaches the eye: half at the mirror on the way in, half on the way
back. Light that the polarizer does not sort perfectly passes on the direct path and forms a ghost.

## Building the eyepiece

We use a plano-convex lens of N-BK7 with a radius of 60 mm, 6 mm thick and 40 mm in diameter. Every
component spawns at the origin facing +y, so the display lies at negative y and the eye at positive y.

````@example pancake_lens
using GLMakie, BeamletOptics, LinearAlgebra, AbstractTrees
GLMakie.activate!(; ssao = true) #hide
const BMO = BeamletOptics
const mm = 1e-3

λ0 = 550e-9
NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)
````

The faces of a lens are its [`FrontFace`](@ref), its [`BackFace`](@ref) and its [`EdgeFace`](@ref).
The films of a [`SurfaceStack`](@ref) are listed from the glass outward: the retarder lies on the
glass, the polarizer on top of it. Their axes are given in the local frame of the lens, so they turn
with it. The retarder is a zero-order plate designed for `λ0`, its retardance `Γ` a quarter wave.

````@example pancake_lens
function pancake_lens(; Γ = π / 2)
    lens = SphericalLens(60mm, Inf, 6mm, 40mm, NBK7)
    surface_model!(lens, FrontFace(), BeamsplitterCoating())
    films = SurfaceStack(Retarder(Γ, [1, 0, -1]; λ0), ReflectivePolarizer([1, 0, 0]))
    surface_model!(lens, BackFace(), films)
    return lens
end
````

The eye looks through the eyepiece from 15 mm behind it. We model it as an ideal lens of 17 mm focal
length and its retina as a [`Detector`](@ref). The display, 10.5 mm in front of the lens near the
focal plane of the folded eyepiece, is only drawn: a [`NonInteractableObject`](@ref) is ignored by the
tracer.

````@example pancake_lens
y_display = -10.5mm
y_eye = 21mm

function eyepiece(; Γ = π / 2)
    lens = pancake_lens(; Γ)
    eye = ThinLens(17mm, 17mm, 12mm, ConstantRefractiveIndex(1.5))
    translate3d!(eye, [0, y_eye, 0])
    retina = Detector(12mm)
    translate3d!(retina, [0, y_eye + 17mm, 0])
    panel = BMO.CuboidMesh(16mm, 1mm, 12mm)
    translate3d!(panel, [-8mm, y_display - 1mm, -6mm])
    return System([lens, eye, retina, NonInteractableObject(panel)])
end

system = eyepiece()
````

## The display

Each pixel emits a cone of circularly polarized rays. A [`PolarizedRay`](@ref) carries its field,
which must be transverse to its direction: for a ray along `d`, the field `(e₁ + i e₂)/√2` with two
transverse unit vectors is circularly polarized.

````@example pancake_lens
function pixel(x, z; λ = λ0, θ_max = deg2rad(10), num_rings = 3, num_azimuths = 8)
    directions = [[0.0, 1.0, 0.0]]
    for k in 1:num_rings, j in 1:num_azimuths
        θ, φ = k / num_rings * θ_max, 2π * j / num_azimuths
        push!(directions, [sin(θ) * cos(φ), cos(θ), sin(θ) * sin(φ)])
    end
    return map(directions) do d
        e₁ = normalize(cross(d, [0, 0, 1]))
        e₂ = cross(e₁, d)
        Beam(PolarizedRay([x, y_display, z], d, λ, (e₁ + im * e₂) / √2))
    end
end
````

Each partial mirror splits a beam into a transmitted and a reflected child, so every cone becomes a
forest of beam trees. `depth_max` bounds the number of splits along a branch.

````@example pancake_lens
function trace!(system, beams; depth_max = 12)
    foreach(beam -> solve_system!(system, beam; depth_max), beams)
    return beams
end

cones = [trace!(system, pixel(x, 0.0; num_rings = 2, num_azimuths = 6)) for x in (-5mm, 0.0, 5mm)];
nothing #hide
````

## Sorting the branches

A branch is a leaf of a beam tree and its ancestors; its rays form one path of the light. The tree
follows the interface of [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl),
whose `Leaves` iterates over the leaves. We count the passes of a branch through the pancake lens: the
image crosses it three times, the ghost once. Branches that do not end on the retina are lost, e.g.
towards the display.

````@example pancake_lens
function branch_rays(leaf)
    path = []
    beam = leaf
    while !isnothing(beam)
        pushfirst!(path, rays(beam))
        beam = AbstractTrees.parent(beam)
    end
    return reduce(vcat, path)
end

function lens_passes(leaf)
    rs = branch_rays(leaf)
    n_lens = NBK7(BMO.wavelength(first(rs)))
    return count(r -> BMO.refractive_index(r) ≈ n_lens, rs)
end

function on_retina(leaf)
    hit = BMO.intersection(last(rays(leaf)))
    return !isnothing(hit) && BMO.object(hit) isa Detector
end

image_branches(beam) = [l for l in Leaves(beam) if on_retina(l) && lens_passes(l) == 3]
ghost_branches(beam) = [l for l in Leaves(beam) if on_retina(l) && lens_passes(l) == 1]

flux(ray) = BMO.refractive_index(ray) * norm(BMO.polarization(ray))^2
````

## The folded path

[`render!`](@ref) draws the components and any ray. We draw the rays of the image branches of three
pixels in a row, each in its own color. Inside the lens the first pass and the reflection at the polarizer overlap with
the second pass; behind the lens each pixel leaves as a collimated beam, which the eye focuses onto
one point of the retina.

````@example pancake_lens
colors = Makie.wong_colors()

function render_branches!(ax, branches; color)
    for leaf in branches, ray in branch_rays(leaf)
        render!(ax, ray; color)
    end
end

fig = Figure(size = (900, 600))
ax = LScene(fig[1, 1]; show_axis = false)
render!(ax, system)
for (i, cone) in enumerate(cones), beam in cone
    render_branches!(ax, image_branches(beam); color = colors[mod1(i, 7)])
end
set_view(ax, [70mm, -30mm, 75mm], [0, 12mm, 0], [0, 0, 1]) #hide
fig
````

## Throughput

The power on the retina in the image, relative to the power of the display:

````@example pancake_lens
display_power = sum(flux(first(rays(beam))) for cone in cones for beam in cone)
image_power = sum(flux(last(rays(l))) for cone in cones for beam in cone for l in image_branches(beam))
image_power / display_power
````

Close to the ideal quarter. At oblique incidence the films see the projections of their axes, which
lets a little light take the direct path.

## The ghost image

A real retarder differs from a quarter wave. With the retardance error δ, a fraction sin²(δ/2) of the
light that enters the lens passes the polarizer on the direct path. The ghost crosses the lens only
once, without the power of the folded mirror, so it does not leave collimated. We trace the on-axis
pixel through an eyepiece whose retarder is off by 25° and draw the image in blue and the ghost in red.

````@example pancake_lens
δ = deg2rad(25)
system_δ = eyepiece(; Γ = π / 2 + δ)
cone_δ = trace!(system_δ, pixel(0.0, 0.0))

fig = Figure(size = (900, 500))
ax = LScene(fig[1, 1]; show_axis = false)
render!(ax, system_δ)
for beam in cone_δ
    render_branches!(ax, image_branches(beam); color = :dodgerblue)
    render_branches!(ax, ghost_branches(beam); color = :red)
end
set_view(ax, [90mm, -10mm, 25mm], [0, 14mm, 0], [0, 0, 1]) #hide
fig
````

On the retina the image of a pixel is a point, the ghost a blurred halo around it. The end of a branch
on the retina, at `y_eye + 17mm`, is the hit point of its last ray:

````@example pancake_lens
function retina_spots(branches)
    points = [BMO.hit_point(last(rays(l))) for l in branches]
    return [p[1] / mm for p in points], [p[3] / mm for p in points]
end

fig = Figure(size = (520, 480))
ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x [mm]", ylabel = "z [mm]",
    title = "What the eye sees: image and ghost of one pixel")
dense_cone = trace!(system_δ, pixel(0.0, 0.0; num_rings = 8, num_azimuths = 32))
for (branches, color, label) in ((ghost_branches, :red, "ghost"), (image_branches, :dodgerblue, "image"))
    x, z = retina_spots(reduce(vcat, branches.(dense_cone)))
    scatter!(ax, x, z; color, markersize = 4, label)
end
axislegend(ax)
fig
````

## Ghost versus retardance error

On axis, the ghost carries `sin²(δ/2) / 2` of the light. The image carries `cos²(δ/2) cos²(δ) / 4`:
the polarizer reflects the fraction `cos²(δ/2)` on the first pass, and the retarder, passed twice more,
turns the fraction `cos²(δ)` into the polarization that the polarizer transmits. Hence
`ghost / image = 2 tan²(δ/2) / cos²(δ)`. A single on-axis ray per eyepiece is enough to check this:

````@example pancake_lens
function ghost_ratio(Γ; λ = λ0)
    beams = trace!(eyepiece(; Γ), pixel(0.0, 0.0; λ, num_rings = 0))
    image = sum(flux(last(rays(l))) for b in beams for l in image_branches(b); init = 0.0)
    ghost = sum(flux(last(rays(l))) for b in beams for l in ghost_branches(b); init = 0.0)
    return ghost / image
end

δs = range(0, deg2rad(30); length = 13)

fig = Figure(size = (600, 380))
ax = Axis(fig[1, 1]; xlabel = "retardance error δ [°]", ylabel = "ghost / image",
    title = "Ghost of a pancake lens")
lines!(ax, rad2deg.(δs), @. 2tan(δs / 2)^2 / cos(δs)^2; label = "analytic", color = :black)
scatter!(ax, rad2deg.(δs), [ghost_ratio(π / 2 + δ) for δ in δs]; label = "traced")
axislegend(ax; position = :lt)
fig
````

## Color

A zero-order quarter-wave plate is a quarter wave at its design wavelength only: at the wavelength λ
its retardance is `π/2 λ0/λ`. Blue and red light therefore leak into the ghost, which is why pancake
lenses use achromatic retarders made of several stacked films.

````@example pancake_lens
λs = range(450e-9, 650e-9; length = 21)

fig = Figure(size = (600, 380))
ax = Axis(fig[1, 1]; xlabel = "wavelength [nm]", ylabel = "ghost / image",
    title = "Ghost of a zero-order retarder designed for 550 nm")
lines!(ax, λs / 1e-9, [ghost_ratio(π / 2; λ) for λ in λs]; color = :black)
scatter!(ax, λs / 1e-9, [ghost_ratio(π / 2; λ) for λ in λs]; color = [wavelength_color(λ) for λ in λs])
fig
````

## Where to go from here

- Replace the single retarder by an achromatic stack of retarders in a [`SurfaceStack`](@ref) and
  compare the color of the ghost.
- The films act wherever they sit: move them to other faces with [`surface_model!`](@ref), or coat
  only a part of a face with a [`MaskedFace`](@ref).
- The [Coatings](@ref) page lists the coating models, including thin-film designs.

