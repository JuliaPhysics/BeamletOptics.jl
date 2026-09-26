```@setup live_michelson
dir = joinpath(@__DIR__, "..", "assets", "examples")

Main.DocUtils.conditional_include(joinpath(dir, "live_michelson_showcase.jl"))
```

# Interactive Michelson interferometer

This example shows how to build an interactive application with [`live_view`](@ref), which is based on the [live rendering](@ref "Live rendering") functions of this package. The Michelson interferometer of the [Michelson interferometer](@ref) tutorial is rendered into a `GLMakie` window, in which the components can be moved and rotated with the mouse and keyboard. After each change, the system is solved again and the beam path, the fringe pattern on the detector and the optical power are updated live.

![Interactive Michelson interferometer](live_michelson_fringes.png)

The full script can be found [here](https://github.com/JuliaPhysics/BeamletOptics.jl/blob/master/docs/src/assets/examples/live_michelson.jl). Running it from the REPL opens the interactive window.

## Setting up the system

The system is identical to the tutorial, but without the optomechanical parts. A [`GaussianBeamlet`](@ref) is used as the HeNe laser source.

```julia
using GLMakie, BeamletOptics

const BMO = BeamletOptics
const cm = 1e-2
const mm = 1e-3

λ = 632.8e-9
w0 = 0.65mm / 2
M2 = 0.7e-3 / BMO.divergence_angle(λ, w0, 1)
beam = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0; M2)

NBK7 = DiscreteRefractiveIndex([λ], [1.51509])

rpm = RightAnglePrismMirror(25mm, 25mm)
zrotate3d!(rpm, deg2rad(45))
translate3d!(rpm, [0, 23.5cm, 0])

cbs = CubeBeamsplitter(BMO.inch, NBK7)
zrotate3d!(cbs, deg2rad(-90))
translate_to3d!(cbs, [18.81cm, 23.5cm, 0])

m1 = RoundPlanoMirror(BMO.inch, 5mm)
zrotate3d!(m1, deg2rad(-90))
translate3d!(m1, [42.715cm, 23.5cm, 0])

m2 = RoundPlanoMirror(BMO.inch, 5mm)
translate3d!(m2, [18.81cm, 37.405cm, 0])

pd_size = 5mm
pd = Detector(pd_size)
translate_to3d!(pd, [18.81cm, 9.595cm, 0])

system = System([rpm, cbs, m1, m2, pd])
```

## Custom updates

The `on_change` callback is called after each solve with the moved component, or `nothing` for the initial solve. Here, it records the optical power on the detector. Since the callback already runs for the initial solve inside [`live_view`](@ref), it only updates an `Observable`, which is plotted after the window has been created. The intensity is evaluated on the full detector area, which makes the movement of the fringes visible:

```julia
full_area = (; x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
power = Observable(Point2f[])

function record_power!(gui, obj)
    P = isnothing(BMO.hits(pd)) ? 0.0 : optical_power(pd; n = 100, full_area...)
    n = isempty(power[]) ? 1 : last(power[])[1] + 1
    push!(power[], Point2f(n, 1e3 * P))
    length(power[]) > 300 && popfirst!(power[])
    notify(power)
    return nothing
end
```

Errors in the callback are logged once and do not interrupt the interaction.

## Opening the interactive window

A single call of [`live_view`](@ref) opens a complete interactive window for the system and the beam: the 3D view, one panel per [`Detector`](@ref) and a status line. The detector panel shows the intensity for Gaussian beamlets and the spot diagram for rays, together with the optical power or the number of rays in its title. By default, the intensity is cropped around the beam, here the full detector area is evaluated instead. The optical power is plotted into an additional axis below the detector panel:

```julia
gui = live_view(system, beam; size = (1200, 700), detectors = [pd => (:intensity, full_area)],
    on_change = record_power!)
power_ax = Axis(gui.fig[1, 2][2, 1]; title = "Optical power", xlabel = "Update", ylabel = "P [mW]")
lines!(power_ax, power; color = :red)
on(_ -> autolimits!(power_ax), power)
display(gui)
```

After each change, `live_view` empties all detectors, solves the system again and updates the beam and the detector panels. There is no need to call [`solve_system!`](@ref) or [`update_render!`](@ref) manually. The figure, the 3D view and the controls are available as `gui.fig`, `gui.ax` and `gui.controls`, e.g. to add static context via `render!(gui.ax, ...)`. Several systems can be shown in the same view via `live_view(system1 => beam1, system2 => beam2)`, and sliders for custom parameters can be added via the `sliders` keyword argument.

## Controls

The components are moved via [`kinematic_controls!`](@ref), see [Interactive kinematics](@ref) for all controls. A click selects a component, a drag on the selected component moves it in the horizontal plane, and every other drag rotates the camera. The arrow keys and `Page Up`/`Page Down` move the selected component along the green, red and blue arrow above it, or rotate it around the rings after switching to the rotate mode with `m`. The step size is changed with `+` and `-`, `Backspace` resets the component and `h` shows an overlay of all controls. The laser can be moved and tilted as well, via the orange marker at its start point. The key `v` switches to the spectator mode, in which the system can be viewed without moving anything by accident. Keyword arguments such as the initial `fine_step` or the `rotation_axis` are passed from `live_view` to [`kinematic_controls!`](@ref).

Since the interferometer is sensitive to changes in the order of the wavelength, the keyboard controls are best suited for alignment. Rotating the mirror `m1` by 1 mrad generates the fringes shown above. Moving the mirror `m2` by ``\lambda/2`` changes the optical path length by ``\lambda``, which corresponds to one period of the optical power:

![Optical power](live_michelson_power.png)
