# Own cards, widgets and controls in the live view

Recipes for `gui = live_view(...)` (needs GLMakie and a window; not for headless scripts). Use them
when adding a component type whose parameters should be editable in the window, or a GUI extension.
Everything works identically in `layout = :compact` and `layout = :app`. Look up any signature with
`scripts/api_lookup.jl <Name>` before use.

## Which path?

| I want to add ... | path | function |
|:--|:--|:--|
| widgets for something in the scene: a component, group, source, system, own `AbstractObject` subtype | a card | `card_rows`, `card_actions` |
| a parameter without a scene object (a setting of the whole setup) | a controls section | `add_controls!` |
| an action (button, or toggle, optionally with a key) | a tool | `add_tool!` |
| a plot or result updated after each solve | a panel | `add_panel!` (see `VISUALIZATION.md`) |

A card needs a scene object: it gives the card its anchor (next to the bounding box, follows the
object), its selection, its pin and its row in the object tree. Without such an object use controls.

## Recipe: a card for your type

The card of the selected object shows the rows of `card_rows` below the head with the buttons of
`card_actions`; only the card of the selection additionally has the step, the Move/Rotate control
and a collapsed "Properties" part (added by the live view, nothing to do).

1. Define `BeamletOptics.card_rows(x::MyType)`. For movable objects start with `pose_card_rows(x)...`.
2. One `CardRow` per line; cells are strings or `CardWidget(T; name, value, on, solve, attributes...)`
   with `T` a Makie block (`Label`, `Slider`, `Toggle`, `Textbox`, `Button`, `Menu`) or an own widget type.
3. `value(gui, obj)` returns what the widget shows, in display units (mm, mrad, %); it must be cheap,
   it runs after every move, solve and input.
4. `on(gui, obj, v)` applies the input `v` to `obj` in SI units. `solve = true` if the input changes
   the optics (a running solve is cancelled first, the systems are solved again afterwards).
5. `name = :something` finds a widget in tests. Extra buttons in the head: a `card_actions(x)` method
   that extends the default ones via `invoke`, as below.
6. No state in the widgets and no captured `gui`: the card may be rebuilt at any time and shows other
   objects with the same declarations.
7. Test it: `gui = live_view(System([...]) => beam; layout)`, `gui.controls.selected[] = obj`, find the
   widget by `name` on the card, set an input, check the object, the shown value and the solve. The
   examples of this file run as tests in `test/Rendering/TestLiveWidgetRecipe.jl`.

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

## Recipe: a card for your system type

A system (subtype of `AbstractSystem`) has a card too: shown, without gizmo and without a selection
(`gui.controls.selected[]` stays `nothing`), by its entry in the component menu (compact) or its row in
the tree (app); `Esc`, a click on empty space or selecting an object closes it. Default rows: number
of objects, number of rays, duration of the last solve. Same steps as above, but no pose. Extend the
default via `invoke`. Objects that are not movable are shown the same way (pose inputs rejected).

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

## Recipe: your own widget type

Any Makie block, or any type constructed as `T(position; attributes...)` that places itself at a
`GridPosition`, can be a `CardWidget` and be used in `add_controls!`, with two methods:

- `card_input(w)` returns the `Observable` whose updates are the user's inputs (each is passed as `v`
  to `on(gui, obj, v)`), or `nothing` (default) for a widget without input.
- `card_show!(w, v)` shows the value `v` of `value(gui, obj)`. Showing is not an input: the card
  ignores the updates of `card_input` while it shows values, but nothing does outside of a card. So do
  not change the observable of `card_input` in `card_show!`, and keep what the user types (the
  `Textbox` method leaves a focused box alone unless `force = true`).

A widget of several blocks must delete them in `Base.delete!` (`@forwarded_layout` does it).
`Makie.@Block` outside of Makie needs `using Makie: make_block_docstring`.

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

1. `add_controls!(gui, title) do layout ... end` builds Makie blocks into `layout`; the layout places the
   section (compact: row above the status row; app: left sidebar). Never use fixed `gui.fig[...]`
   positions.
2. A `Textbox` or `Menu` built inside `f` takes the keyboard automatically (3D keys ignored while
   focused or open); blocks added after `f` returns do not.
3. Every change of the optics goes through `retrace!(f, gui)`, which cancels a running solve first and
   then calls `f()`. Never change objects outside `f`.
4. Return `update(gui)` from the builder to show the current state after each full solve and once right
   away (e.g. after an undo). `update` must not change objects. A builder ending with any other value,
   including the `on(...)` of a listener, has no `update`.
5. Errors in callbacks are logged once and do not stop the view.
6. Identical in `:compact` and `:app`; own widget types use `card_input`/`card_show!` here too.

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

## Pitfalls

- Do not put a parameter without a scene object on a card; there is no `add_controls!(gui, title, obj)`.
- `value` in display units, `on` in SI units; forgetting the conversion (mm to m, mrad to rad) is the
  usual bug.
- `on` that changes the optics needs `solve = true`; without it nothing is solved.
- Changes of objects from a controls widget or tool go through `retrace!(f, gui)`, not directly.
- `card_show!` must not change the observable of `card_input` (it would count as an input).
- A `do` block of `add_controls!` or `add_panel!` that ends with `on(...)` has no `update`; end it with
  `return gui -> ...` when one is wanted.
- A card keeps no state: do not capture `gui` or the object in the widget declarations.

## Checklist for a new component

- `interact3d(system, object, beam, ray)`: required.
- `render!`: optional, how it is drawn.
- `card_rows` (and `card_actions`): optional, its rows and buttons on the card.
- `BeamletOptics.properties`: optional, the properties in the "Properties" part of the card.
