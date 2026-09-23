using GLMakie, BeamletOptics

const BMO = BeamletOptics
const cm = 1e-2
const mm = 1e-3

## Michelson interferometer, see the Michelson tutorial
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

## Optical power over the updates, plotted into an additional axis below the detector panel
power = Observable(Point2f[])
power_ax = nothing

"""Returns the optical power [W] from the intensity of a detector panel of the `live_view`."""
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

## Interactive window, the detector intensity is evaluated on the full detector area
full_area = (; x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
gui = live_view(system, beam; size = (1200, 700), detectors = [pd => (:intensity, full_area)],
    on_change = record_power!)
fig = gui.fig
controls = gui.controls

# Solves the system again after moving a component from code, as the controls do after each change
on_change(obj) = controls.on_change(obj)

# Open the interactive window when used from the REPL or run as a script
if isinteractive()
    display(gui)
elseif abspath(PROGRAM_FILE) == @__FILE__
    wait(display(gui))
end
