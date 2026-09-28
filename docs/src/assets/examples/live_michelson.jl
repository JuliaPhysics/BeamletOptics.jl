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

## Optical power over the updates, the detector intensity is evaluated on the full detector area
full_area = (; x_min = -pd_size / 2, x_max = pd_size / 2, z_min = -pd_size / 2, z_max = pd_size / 2)
power = Point2f[]

# Called after each full solve, also while the power panel is hidden (a tab of the app layout)
function record_power!(gui, obj)
    P = isnothing(BMO.hits(pd)) ? 0.0 : optical_power(pd; n = 100, full_area...)
    n = isempty(power) ? 1 : last(power)[1] + 1
    push!(power, Point2f(n, 1e3 * P))
    length(power) > 300 && popfirst!(power)
    return nothing
end

## Interactive window, `layout = :app` opens it as an application window
gui = live_view(system, beam; size = (1200, 700), detectors = [pd => (:intensity, full_area)],
    on_change = record_power!, layout = :compact)

# The optical power as an own panel: below the detector panel, or a tab in the app layout
add_panel!(gui, "Optical power") do layout
    ax = Axis(layout[1, 1]; xlabel = "Update", ylabel = "P [mW]")
    pts = Observable(copy(power))
    lines!(ax, pts; color = :red)
    # Called after each full solve while the panel is shown
    return gui -> (pts[] = copy(power); autolimits!(ax))
end
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
