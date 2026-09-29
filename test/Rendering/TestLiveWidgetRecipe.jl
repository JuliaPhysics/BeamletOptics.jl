module TestLiveWidgetRecipe

#=
The code examples of the recipe for the widgets of the live view (the docs page "Live view widgets"
and `skills/beamletoptics/WIDGETS.md`), run as tests. Each recipe is a block of definitions between
two `# ---` lines, which the docs copy unchanged, followed by the tests of what the recipe promises.
=#

using BeamletOptics
using Makie
using Test
# exported with the release of the recipe; imported explicitly until then
using BeamletOptics: card_input, card_show!

#=
Shared by the tests: a beam along +y, a mirror at 45° that reflects it along +x onto a detector
=#

function _bench()
    mirror = RoundPlanoMirror(25e-3, 5e-3)
    zrotate3d!(mirror, deg2rad(45))
    translate3d!(mirror, [0, 0.1, 0])
    pd = Detector(5e-3)
    zrotate3d!(pd, -π / 2)
    translate3d!(pd, [0.1, 0.1, 0])
    return mirror, pd
end

_hits(pd) = something(BeamletOptics.hits(pd), []) |> length

# --- Recipe: a card for your type ---

# A filter that passes the rays and attenuates their power by `transmission` (0-1)
mutable struct MyAttenuator{T} <: BeamletOptics.AbstractObject{T}
    shape::BeamletOptics.PlanoSurfaceSDF{T}
    transmission::Float64
end

MyAttenuator(diameter, thickness; transmission = 0.5) =
    MyAttenuator(BeamletOptics.PlanoSurfaceSDF(thickness, diameter), Float64(transmission))

# The rays pass straight through (the power is not traced in this example)
function BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, a::MyAttenuator,
        ::Beam{T, R}, ray::R) where {T <: Real, R <: Ray{T}}
    pos = position(ray) + length(ray) * direction(ray)
    # entering: the next surface is the back of the filter
    hint = BeamletOptics.isentering(ray) ? BeamletOptics.Hint(a) : nothing
    new = Ray{T}(pos, direction(ray), nothing, BeamletOptics.wavelength(ray), BeamletOptics.refractive_index(ray))
    return BeamletOptics.BeamInteraction{T, R}(hint, new)
end

# The card: the pose, a slider for the transmission with its value in %, and an action "block"
# next to the default ones ("hide")
BeamletOptics.card_rows(a::MyAttenuator) = (pose_card_rows(a)...,
    CardRow("T",
        CardWidget(Slider; name = :transmission, range = 0:0.01:1, width = 150,
            value = (gui, a) -> a.transmission,
            on = (gui, a, v) -> (a.transmission = v), solve = true),
        CardWidget(Label; name = :percent,
            value = (gui, a) -> "$(round(Int, 100 * a.transmission)) %")))

BeamletOptics.card_actions(a::MyAttenuator) = (
    invoke(BeamletOptics.card_actions, Tuple{BeamletOptics.AbstractObject}, a)...,
    CardWidget(Button; name = :block, label = "block",
        on = (gui, a, _) -> (a.transmission = 0.0), solve = true))

# ---

# --- Recipe: a card for your system type ---

# A system with a name, traced like a `System`
struct MyBench <: BeamletOptics.AbstractSystem
    objects::Vector{BeamletOptics.AbstractObject}
    name::String
end

BeamletOptics.objects(b::MyBench) = b.objects

# The default rows of a system (objects, rays, solve) and the name of the bench
BeamletOptics.card_rows(b::MyBench) = (
    invoke(BeamletOptics.card_rows, Tuple{BeamletOptics.AbstractSystem}, b)...,
    CardRow("bench", CardWidget(Label; name = :bench, value = (gui, b) -> b.name)))

# ---

# --- Recipe: your own widget type ---

# `Makie.@Block` outside of Makie needs this name of Makie
using Makie: make_block_docstring

# A value with a "−" and a "+" button, which step it by `step`
Makie.@Block Stepper begin
    @forwarded_layout
    minus::Button
    plus::Button
    @attributes begin
        "The value shown, see `card_show!`."
        value::Float64 = 0.0
        "The value requested by the last click, see `card_input`."
        input::Float64 = 0.0
        "The change of the value per click."
        step = 1.0
        "The horizontal alignment of the block in its suggested bounding box."
        halign = :center
        "The vertical alignment of the block in its suggested bounding box."
        valign = :center
        "The width setting of the block."
        width = Auto()
        "The height setting of the block."
        height = Auto()
        "Controls if the parent layout can adjust to this block's width."
        tellwidth::Bool = true
        "Controls if the parent layout can adjust to this block's height."
        tellheight::Bool = true
        "The align mode of the block in its parent GridLayout."
        alignmode = Inside()
    end
end

function Makie.initialize_block!(s::Stepper)
    s.minus = Button(s.layout[1, 1]; label = "−")
    Label(s.layout[1, 2], lift(v -> string(round(v; digits = 3)), s.value); width = 50)
    s.plus = Button(s.layout[1, 3]; label = "+")
    # a click is an input: the value one step away from the shown one
    on(_ -> (s.input[] = s.value[] - s.step[]), s.minus.clicks)
    on(_ -> (s.input[] = s.value[] + s.step[]), s.plus.clicks)
    return nothing
end

# The protocol of the widgets of the cards: the inputs, and showing a value (not an input)
BeamletOptics.card_input(s::Stepper) = s.input
BeamletOptics.card_show!(s::Stepper, v) = (s.value[] = v; nothing)

# A beam block, which its card moves up and down (along z) in steps
struct MyBeamBlock{T} <: BeamletOptics.AbstractObject{T}
    shape::BeamletOptics.PlanoSurfaceSDF{T}
end

MyBeamBlock(width, thickness) = MyBeamBlock(BeamletOptics.PlanoSurfaceSDF(thickness, width))

# Absorbs the rays that hit it
BeamletOptics.interact3d(::BeamletOptics.AbstractSystem, ::MyBeamBlock, ::Beam, ::Ray) = nothing

BeamletOptics.card_rows(b::MyBeamBlock) = (CardRow("z [mm]",
    CardWidget(Stepper; name = :height, step = 15.0,
        value = (gui, b) -> 1e3 * position(b)[3],
        on = (gui, b, z) -> translate_to3d!(b, [position(b)[1], position(b)[2], 1e-3 * z]),
        solve = true)),)

# The same widget in a controls section, for the block: shows its position after each solve
function add_block_controls!(gui, block)
    add_controls!(gui, "Beam block") do layout
        Label(layout[1, 1], "z [mm]")
        stepper = Stepper(layout[1, 2]; step = 15.0)
        on(card_input(stepper)) do z
            retrace!(() -> translate_to3d!(block, [position(block)[1], position(block)[2], 1e-3 * z]), gui)
        end
        return gui -> card_show!(stepper, 1e3 * position(block)[3])
    end
end

# ---

# --- Recipe: controls without a scene object ---

# A misalignment scenario of the setup, not a property of an object: it tilts the fold mirror about z;
# the tilt [rad] and whether it is applied
const tilt = Ref(0.0)
const tilted = Ref(false)

function set_tilt!(gui, mirror, θ, on::Bool)
    retrace!(gui) do
        # the rotation from the applied tilt to the new one
        zrotate3d!(mirror, (on ? θ : 0.0) - (tilted[] ? tilt[] : 0.0))
        tilt[], tilted[] = θ, on
    end
    return nothing
end

function add_tilt_controls!(gui, mirror)
    add_controls!(gui, "Alignment") do layout
        Label(layout[1, 1], "tilt [mrad]")
        box = Textbox(layout[1, 2]; placeholder = "0", width = 80, validator = Float64)
        toggle = Toggle(layout[1, 3])
        Label(layout[1, 4], "tilted")
        on(s -> set_tilt!(gui, mirror, 1e-3 * parse(Float64, s), toggle.active[]), box.stored_string)
        on(v -> v == tilted[] || set_tilt!(gui, mirror, tilt[], v), toggle.active)
        # the parameter as it is, e.g. after it was set from code
        return gui -> (card_show!(box, 1e3 * tilt[]); card_show!(toggle, tilted[]))
    end
end

# ---

#=
Tests
=#

const Ext = Base.get_extension(BeamletOptics, :BeamletOpticsMakieExt)

# A live view whose changes are solved at once; `on_change` counts the full solves
function _live_view(pairs...; kwargs...)
    solves = Ref(0)
    gui = live_view(pairs...; trace_budget = Inf, throttle = false,
        on_change = (gui, obj) -> (solves[] += 1), kwargs...)
    return gui, solves
end

# The widget `name` on the card of the object of the card of the selection: floating (compact) or
# docked in the inspector (app)
_card(gui) = gui.cards.selection
_card(gui::Ext.LiveView{Ext.AppLayout}) = gui.layout.inspector.card
_widget(gui, name) = Ext._card_widget(_card(gui), name)

@testset "Live view widget recipes" begin
    @test !isnothing(Ext)

    @testset "a card for your type ($layout)" for layout in (:compact, :app)
        mirror, pd = _bench()
        a = MyAttenuator(25e-3, 2e-3; transmission = 0.3)
        translate3d!(a, [0, 0.05, 0])
        gui, solves = _live_view(System([a, mirror, pd]) => Beam([0.0, 0, 0], [0.0, 1, 0]); layout)
        # the rays pass the attenuator
        @test _hits(pd) == 1
        gui.controls.selected[] = a
        # the pose rows and the declared widgets, found by their names
        @test !isnothing(_widget(gui, :x))
        slider = _widget(gui, :transmission)
        @test slider isa Slider
        @test slider.value[] ≈ 0.3
        @test _widget(gui, :percent).text[] == "30 %"
        # an input is applied and solved
        n = solves[]
        Makie.set_close_to!(slider, 0.8)
        @test a.transmission ≈ 0.8
        @test _widget(gui, :percent).text[] == "80 %"
        @test solves[] == n + 1
        # the extra action next to the default ones
        @test _widget(gui, :hide) isa Button
        notify(_widget(gui, :block).clicks)
        @test a.transmission == 0
        @test slider.value[] == 0
        @test solves[] == n + 2
        close(gui)
    end

    @testset "a card for your system type ($layout)" for layout in (:compact, :app)
        mirror, pd = _bench()
        bench = MyBench([mirror, pd], "Fold")
        gui, _ = _live_view(bench => Beam([0.0, 0, 0], [0.0, 1, 0]); layout)
        @test _hits(pd) == 1
        # the system is inspected via its entry in the component menu or its row in the tree
        h = only(gui.system_handles)
        layout == :compact ? (gui.widgets.menu.i_selected[] = 1) : (gui.layout.tree.clicked[] = h)
        @test gui.objects.inspected === bench
        @test isnothing(gui.controls.selected[])
        @test _widget(gui, :bench).text[] == "Fold"
        @test _widget(gui, :objects).text[] == "2"
        @test _widget(gui, :rays).text[] == "1"
        close(gui)
    end

    @testset "your own widget type ($layout)" for layout in (:compact, :app)
        mirror, pd = _bench()
        block = MyBeamBlock(20e-3, 2e-3)
        zrotate3d!(block, -π / 2)
        translate3d!(block, [0.05, 0.1, 0])
        gui, solves = _live_view(System([mirror, pd, block]) => Beam([0.0, 0, 0], [0.0, 1, 0]); layout)
        # the block is in the beam
        @test _hits(pd) == 0
        # on a card: the value shown, a click is an input, which is solved
        gui.controls.selected[] = block
        s = _widget(gui, :height)
        @test s isa Stepper
        @test s.value[] ≈ 0
        n = solves[]
        notify(s.plus.clicks)
        @test position(block)[3] ≈ 0.015
        @test s.value[] ≈ 15
        @test solves[] == n + 1
        # above the beam, which hits the detector again
        @test _hits(pd) == 1
        # in a controls section: shows the position right away and after each solve
        layout_ = add_block_controls!(gui, block)
        stepper = only(b for b in Ext._blocks!(Any[], layout_) if b isa Stepper)
        @test stepper.value[] ≈ 15
        notify(stepper.plus.clicks)
        @test position(block)[3] ≈ 0.03
        @test stepper.value[] ≈ 30
        @test s.value[] ≈ 30
        # moved on the card: the controls follow after the solve
        notify(s.minus.clicks)
        @test position(block)[3] ≈ 0.015
        @test stepper.value[] ≈ 15
        # back into the beam, which the solve detects: it traces from the start of the beam
        notify(s.minus.clicks)
        @test position(block)[3] ≈ 0 atol = 1e-12
        @test stepper.value[] ≈ 0 atol = 1e-9
        @test _hits(pd) == 0
        # the card is rebuilt for another object: the blocks of the stepper are deleted with it
        @test !isempty(s.blockscene.children)
        gui.controls.selected[] = mirror
        @test isnothing(_widget(gui, :height))
        @test isempty(s.blockscene.children)
        close(gui)
    end

    @testset "controls without a scene object ($layout)" for layout in (:compact, :app)
        tilt[], tilted[] = 0.0, false
        mirror, pd = _bench()
        gui, solves = _live_view(System([mirror, pd]) => Beam([0.0, 0, 0], [0.0, 1, 0]); layout)
        R0 = copy(orientation(mirror))
        controls = add_tilt_controls!(gui, mirror)
        blocks = Ext._blocks!(Any[], controls)
        box = only(b for b in blocks if b isa Textbox)
        toggle = only(b for b in blocks if b isa Toggle)
        # the textbox takes the keyboard
        @test box in gui.custom.boxes
        # shown right away
        @test box.displayed_string[] == "0.0"
        @test !toggle.active[]
        # the inputs change the parameter via `retrace!`, i.e. solved
        n = solves[]
        box.stored_string[] = "20"
        @test tilt[] ≈ 0.02
        @test orientation(mirror) ≈ R0
        toggle.active[] = true
        @test tilted[]
        @test !(orientation(mirror) ≈ R0)
        @test solves[] == n + 2
        # the reflected beam misses the detector
        @test _hits(pd) == 0
        # changed from code: the controls show the parameter after the solve
        set_tilt!(gui, mirror, 0.0, false)
        @test box.displayed_string[] == "0.0"
        @test !toggle.active[]
        @test orientation(mirror) ≈ R0
        @test _hits(pd) == 1
        close(gui)
    end
end

end # module
