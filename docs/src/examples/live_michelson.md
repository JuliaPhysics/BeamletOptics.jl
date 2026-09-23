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

## Opening the interactive window

A single call of [`live_view`](@ref) opens a complete interactive window for the system and the beam: the 3D view, one panel per [`Detector`](@ref) and a status line. The detector panel shows the intensity for Gaussian beamlets and the spot diagram for rays, together with the optical power or the number of rays in its title. By default, the intensity is cropped around the beam. Here, the full detector area is evaluated instead, which makes the movement of the fringes visible. The `record_power!` callback is defined in the [next section](@ref "Custom updates").

```julia
full_area = (; x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
gui = live_view(system, beam; size = (1200, 700), detectors = [pd => (:intensity, full_area)],
    on_change = record_power!)
display(gui)
```

After each change, `live_view` empties all detectors, solves the system again and updates the beam and the detector panels. There is no need to call [`solve_system!`](@ref) or [`update_render!`](@ref) manually. The figure, the 3D view and the controls are available as `gui.fig`, `gui.ax` and `gui.controls`, e.g. to add static context via `render!(gui.ax, ...)`. Several systems can be shown in the same view via `live_view(system1 => beam1, system2 => beam2)`, and sliders for custom parameters can be added via the `sliders` keyword argument.

## Custom updates

The `on_change` callback is called after each solve with the moved component, or `nothing` for the initial solve. Here, it records the optical power on the detector, which is calculated from the intensity of the detector panel, and plots it into an additional axis below the panel:

```julia
power = Observable(Point2f[])
power_ax = nothing

function panel_power(panel)
    isnothing(BMO.hits(panel.pd)) && return 0.0
    x, z, I = panel.heat_x[], panel.heat_y[], panel.heat_I[] # [mm], [mm], [W/m²]
    return sum(I) * (x[2] - x[1]) * (z[2] - z[1]) * mm^2
end

function record_power!(gui, obj)
    global power_ax
    # Called for the first time after the window has been set up
    if isnothing(power_ax)
        power_ax = Axis(gui.fig[1, 2][2, 1]; title = "Optical power", xlabel = "Update", ylabel = "P [mW]")
        lines!(power_ax, power; color = :red)
    end
    n = isempty(power[]) ? 1 : last(power[])[1] + 1
    push!(power[], Point2f(n, 1e3 * panel_power(gui.panels[1])))
    length(power[]) > 300 && popfirst!(power[])
    notify(power)
    autolimits!(power_ax)
    return nothing
end
```

Errors in the callback are logged once and do not interrupt the interaction.

## Controls

The components are moved via [`kinematic_controls!`](@ref), which enables the following controls.
A click selects a component, a drag on the selected component moves or rotates it, and every
other drag rotates the camera as usual, so rotating the camera never selects or moves a component
by accident:

| Input                          | Move mode                    | Rotate mode                  |
|:-------------------------------|:-----------------------------|:-----------------------------|
| Left-click on a component      | Select it                    | Select it                    |
| Left-drag on the selection     | Move in the horizontal plane | Rotate around the blue axis  |
| Left-drag elsewhere            | Rotate the camera            | Rotate the camera            |
| `↑` / `↓`                      | Move along the green arrow   | Rotate around the red ring   |
| `→` / `←`                      | Move along the red arrow     | Rotate around the blue ring  |
| `Page Up` / `Page Down`        | Move along the blue arrow    | Rotate around the green ring |
| `+` / `-`                      | Increase / decrease the step | Increase / decrease the step |

Further controls: `m` switches between the move and the rotate mode, pressing shift multiplies
the step size by 10, `Backspace` resets the selected component to its initial pose, `Esc` or a
click on empty space deselects it and `h` shows or hides an overlay of all controls. The keys `+` and `-`
change the step size of the current mode along the 1-2-5 sequence, e.g. 10 nm → 20 nm → 50 nm →
100 nm, also without a selected component. The current step size is shown in the hint line at the
top of the 3D view. Keyword arguments such as the initial `fine_step` or the `rotation_axis` are
passed from `live_view` to [`kinematic_controls!`](@ref).

The selected component is marked by a box and three axes above it: its local y-axis (green), its
local x-axis (red) and the vertical rotation axis (blue). In the move mode the axes are shown as
arrows, in the rotate mode as rings. The first key of each pair moves the component in the
direction of the arrow, or rotates it in the direction of the ring.

Since the interferometer is sensitive to changes in the order of the wavelength, the keyboard controls are best suited for alignment. Rotating the mirror `m1` by 1 mrad generates the fringes shown above. Moving the mirror `m2` by ``\lambda/2`` changes the optical path length by ``\lambda``, which corresponds to one period of the optical power:

![Optical power](live_michelson_power.png)
