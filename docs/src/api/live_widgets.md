# Live view widgets

This page is for package authors and developers who add a component or a GUI extension to
[`live_view`](@ref). It shows how a new type gets its own card in the 3D view (and in the
"Properties" sidebar of `layout = :app`), how an own widget type is used on cards and in controls,
and how parameters without a scene object get their own controls. Everything works identically in
`layout = :compact` and `layout = :app`. The user-facing description of the window is on the page
[Visualization](@ref).

All code blocks need a `Makie` backend with a window (`GLMakie`) and are therefore not run when the
docs are built. They are run as tests in `test/Rendering/TestLiveWidgetRecipe.jl`, from which they are
copied unchanged.

## Which path?

| I want to add ... | path | function |
|:--|:--|:--|
| widgets for something in the scene: a component, a group, a source, a system, an own `AbstractObject` subtype | a card | [`card_rows`](@ref), [`card_actions`](@ref) |
| a parameter without a scene object, e.g. a setting of the whole setup | a controls section | [`add_controls!`](@ref) |
| an action, e.g. a button or a toggle with a key | a tool | [`add_tool!`](@ref) |
| a plot or a result that is updated after each solve | a panel | [`add_panel!`](@ref) |

A card needs a scene object because the object gives the card its anchor (the card is placed next
to the bounding box of the object and follows it), its selection (the card shows what is selected,
clicked in the menu or in the tree), its pin (the card stays with its object) and its row in the
object tree. A parameter without such an object has none of these, hence it belongs into a controls
section.

## Recipe: a card for your type

The card of a selected object shows the rows of [`card_rows`](@ref) below its head with the actions
of [`card_actions`](@ref), and, only on the card of the selection, the step, the Move/Rotate control
and the collapsed "Properties" part, which the live view adds. To add rows for an own type:

1. Define a method `BeamletOptics.card_rows(x::MyType)`. For a movable object, start from
   `pose_card_rows(x)...` ([`pose_card_rows`](@ref)), which gives the position and rotation boxes.
2. Return one [`CardRow`](@ref) per line. Its cells are strings (labels) or [`CardWidget`](@ref)s.
3. `value(gui, obj)` returns what the widget shows, in display units (mm, mrad, %), and must be
   cheap: it runs after every move, solve and input.
4. `on(gui, obj, v)` applies an input `v` to `obj` in SI units. Set `solve = true` if the input
   changes the optics: a running solve is cancelled first and the systems are solved again after
   the input.
5. Give a widget a `name` to find it in tests. Extra actions in the head of the card are added by a
   method of [`card_actions`](@ref), which extends the default ones ("hide").
6. Keep no state in the widgets and capture no `gui`: the card may be rebuilt at any time and
   shows other objects with the same declarations. `value` and `on` get the `gui` and the object.
7. Test it: create a `live_view` with the object, select it with `gui.controls.selected[] = obj`,
   look up the widget by its `name` on the card, set an input and check the object, the shown value
   and that a solve followed. The tests of the recipes do this for both layouts, see
   `test/Rendering/TestLiveWidgetRecipe.jl`.

The example is a filter with a transmission, whose card has a slider (its value is shown in %) and an
action "block" next to "hide":

```julia
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
```

The types of BeamletOptics have such methods as well, e.g. the ray count slider of a source or the
mode of the panel of a `Detector`.

## Recipe: a card for your system type

A system (a subtype of `AbstractSystem`) has a card as well: it is shown, without a gizmo and
without a selection (`gui.controls.selected[]` stays `nothing`), when its entry in the component menu
(`layout = :compact`) or its row in the object tree (`layout = :app`) is clicked, and closed by
`Esc`, a click on empty space or the selection of an object. By default, it shows the number of its
objects, the number of rays and the duration of the last solve. The steps are those of the card of a
type, except that a system has no pose. Extend the default rows via `invoke`:

```julia
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
```

An object that is not movable is shown on its card in the same way (inspected instead of selected),
its pose boxes reject inputs with a message in the status line.

## Recipe: your own widget type

Any `Makie` block can be a [`CardWidget`](@ref) and is then placed, hidden and refreshed by the card.
An own widget type, e.g. a `Makie.@Block` or another type that is constructed as
`T(position; attributes...)` and places itself at a `GridPosition`, takes part by two methods:

- [`card_input`](@ref)`(w)` returns the `Observable` whose updates are the inputs of the user. Each
  update is passed as `v` to the `on(gui, obj, v)` of the declaration. Return `nothing` (the default)
  if the widget takes no input.
- [`card_show!`](@ref)`(w, v)` shows the value `v` of `value(gui, obj)`. Showing a value is not an
  input: the card ignores the updates of `card_input` while it shows values, but nothing does outside
  of a card, e.g. in [`add_controls!`](@ref). Hence do not change the observable of `card_input`
  in `card_show!`, and keep what the user is typing (the method for a `Textbox` does not change a
  focused box unless `force = true`).

An own widget of several blocks must delete them in `Base.delete!`, which the card calls when it
rebuilds its widgets; a `Makie.@Block` with `@forwarded_layout` does it by itself. `Makie.@Block`
outside of `Makie` needs the name `make_block_docstring` of `Makie` in scope. The example is a value
with a "−" and a "+" button, first on the card of a block, then in a controls section:

```julia
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
```

## Recipe: controls without a scene object

A parameter that belongs to no object, e.g. an alignment setting of the whole setup, has no card.
It is a controls section:

1. `add_controls!(gui, title) do layout ... end` builds `Makie` blocks into the given `layout`. The
   layout places the section (a row above the status row, or a section of the left sidebar); never
   use fixed `gui.fig[...]` positions.
2. A `Textbox` or `Menu` built inside `f` takes the keyboard automatically, i.e. the keys of the 3D
   view are ignored while it is focused or open. Blocks added after `f` returned do not.
3. Every change of the optics goes through `retrace!(f, gui)` ([`retrace!`](@ref)), which cancels a
   running solve first and makes the change `f` in its place. Never change objects outside `f`.
4. Return `update(gui)` from the builder to show the current state after each full solve (and once
   right away), e.g. after an undo or a change from code. `update` must not change objects. A
   builder that ends with any other value, including the `on(...)` of a listener, has no `update`.
5. Errors in callbacks are logged once and do not stop the view.
6. It works identically in `:compact` and `:app`; own widget types use `card_input` and `card_show!`
   here too, as in the example above.

The example is a misalignment scenario, set with a textbox and switched on and off with a toggle. It
acts on the fold mirror, but it is a setting of the setup, not a property of the mirror:

```julia
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
```

## Checklist for a new component

- `interact3d(system, object, beam, ray)`: required, see [Core design](@ref).
- `render!`: optional, how the component is drawn.
- `card_rows` (and `card_actions`): optional, its rows and buttons on the card of the live view.
- `BeamletOptics.properties`: optional, the properties listed in the "Properties" part of the card
  and of the inspector, see [`properties`](@ref).

## Reference

The functions and types of the recipes:

```@docs; canonical=false
card_rows
pose_card_rows
card_actions
CardRow
CardWidget
card_input
card_show!
add_controls!
add_tool!
```
