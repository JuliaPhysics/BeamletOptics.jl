# BeamletOptics.jl (BMO): developer instructions for coding agents

A Julia package for non-sequential 3D ray and Gaussian beamlet tracing, used to simulate
breadboard optical setups (e.g. laser interferometers) with lenses, mirrors, beamsplitters,
polarizers and detectors. See [docs/src/index.md](docs/src/index.md) for the pitch and
[docs/src/api/conventions.md](docs/src/api/conventions.md) for the binding physical and
geometric conventions (global optical axis +y, right-handed frames, CCW rotations).

This file is for **developing** BMO. Guidance for **using** BMO lives in the agent skill
[skills/beamletoptics/](skills/beamletoptics/SKILL.md) (see "Agent skill" below). Code you
write here is read and maintained by humans: follow "Code style" below as strictly as the
physics.

## Design philosophy

BMO is meant to feel like a **digital laboratory**: the user places components in 3D space
the way they would arrange them on an optical breadboard, and the tracer works out the rest.
Evaluate every API and architecture decision against this. The core principles, from
[docs/src/api/core.md](docs/src/api/core.md):

1. Optical interaction is decoupled from geometry representation.
2. Optical elements are closed volumes, or must mimic one (exceptions apply, e.g. coatings).
3. Elements must be freely movable and work for (almost) any angle of incidence: no paraxial
   shortcuts, no assumed canonical orientation.
4. Without extra knowledge, tracing is non-sequential: the solver finds what a ray or beam
   hits next by searching the scene, not from a user-declared path or object order.
5. With extra knowledge (a `Hint`), tracing can go sequential. This is an optimization and a
   tool for component authors, never a requirement placed on the user.

**The extension promise:** a developer defines a new `AbstractObject` subtype and its
`interact3d(system, object, beam, ray)` method (plus `intersect3d` if it needs custom geometry,
optionally `bounding_sphere_of` to make misses of that geometry cheap, and `initialize!` if it
stores data during a solve), and the rest of the API (kinematics, threading) works without
further integration. A component may also add `card_rows` (and `card_actions`) for its card in
the interactive GUI, the separate package
[BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl). When adding
infrastructure, push complexity into the generic solver instead of asking component authors to
handle it. Current exception: coincident boundaries (plate beamsplitters, cemented doublets)
are disambiguated per component by returning a `Hint`
([src/AbstractTypes/AbstractSystem.jl](src/AbstractTypes/AbstractSystem.jl), consumed in
[src/System.jl](src/System.jl)); core.md calls this a burden on the developer.

**Be suspicious of** a component or algorithm that assumes a world-space orientation, an
object order in a `System` or adjacency not derived from geometry at trace time; of paraxial or
near-normal-incidence shortcuts without a documented limit on the angle of incidence; and of
functionality that makes users declare sequence or structure the solver could infer.

## Code map

- `src/AbstractTypes/`: the interfaces (`AbstractObject`, `AbstractShape`, shape traits,
  `AbstractRay`/`Intersection`, `AbstractBeam`, `AbstractSystem`/`Hint`, kinematic trait).
- `src/System.jl`: the intersect-interact loop and tracing.
- `src/Rays.jl`, `PolarizedRays.jl`, `Beam.jl`, `Gaussian.jl`, `AstigmaticGaussian.jl`,
  `BeamGroups/`: beam models and sources.
- `src/OpticalComponents/`: components, one family per folder or file.
- `src/SDFs/`, `src/Mesh.jl`: geometry backends.
- `src/Exports.jl`: the public API (`export` and `public`). Changing it affects the agent skill.
- `ext/`: `BeamletOpticsMakieExt` and its `Render*.jl` files (`render!` and live rendering behind
  the render handle protocol of `src/Render/RenderHandles.jl`; `src/Render/` holds the entry
  points and docstrings, one file per extension file). `BeamletOpticsOpticsBaseExt` couples BMO
  to other solvers through the `PlaneField` of
  [OpticsBase](https://github.com/StackEnjoyer/OpticsBase.jl).
- The GUI uses the exported names, the render handle protocol and the `public` names of
  `src/Exports.jl` (see "Public API for dependent packages" in
  [docs/src/api/api.md](docs/src/api/api.md)): treat them like exported API. It also uses some
  internal names (the abstract types, the kinematic and shape traits, `intersection`, `hits`,
  `objects`): before renaming one, search the GUI for it and update the GUI in step.
- `skills/beamletoptics/`: the user-facing agent skill.

## Workflow for a change

1. **Before writing code**, list what already covers the task: the exported and `public` names
   in `src/Exports.jl`, the functions in the docstring of the abstract type you touch, and the
   verbs of `Base` and AbstractTrees. Name the verb or type you extend (V1). A new public name
   needs a reason: new physics or a new user action.
2. **While writing**, follow "Code style". Stop and ask the maintainers, instead of choosing
   silently, if a requirement can only be met by breaking a rule, if two rules conflict, or if
   the change would alter an interface contract: the `# Implementation reqs.` of an abstract
   type, the meaning of a hook or of a public verb (e.g. "the first beam is the center ray").
3. **Before handing over**, re-read your whole diff once per rule group of "Code style" and fix
   every violation. A change is done when it also has:
   - tests for new or changed behavior (X1);
   - correct docstrings, and docs or examples where users need them (C1);
   - an updated Mermaid diagram if it changes a core system, or a new one for a new core
     system (C1);
   - the agent skill updated if public behavior changed ("Agent skill");
   - the GUI searched for every renamed or removed name.
4. **In the PR text**, state the style check: "Style: checked against AGENTS.md" plus any
   exception you kept and why, and which skill files were updated (or that none were).

## Code style

Guiding idea: **few verbs, many methods, typed dispatch.** A user handles BMO like a lab: take a
part from the catalog, place it, run, measure. New behavior is a new method of an existing verb,
not a new name, flag or mechanism. The [SciML Style Guide](https://github.com/SciML/SciMLStyle)
is the formatting baseline; these rules take precedence.

**Scope:** all rules apply to core code (physics, geometry, solver, components, detectors,
sources, `ext/BeamletOpticsOpticsBaseExt.jl`). Rendering and GUI support code (`src/Render/`,
`ext/Render*.jl`, `ext/BeamletOpticsMakieExt.jl`, `src/Utils/Progress.jl`, `src/Properties.jl`)
follows only V1, V2 and D1 for exported or `public` names.

### V · Verbs and API

**V1 · Extend an existing verb before inventing a name.** Add a method to the verb that already
means this (in BMO, `Base`, AbstractTrees); never rewrite an existing verb or getter inline
(`dot(a, b) / (norm(a) * norm(b))` where `cos_angle3d` exists, a field read where a getter
exists). A different quantity gets its own verb, the same quantity on another level a method of
the same verb; no Bool flag that switches what a function returns. Normalize arguments with a
method that converts and calls the verb again; an optional argument that selects a path is a
positional method on `::Nothing`.
```julia
# Don't
hit_count(d::Detector) = isnothing(hits(d)) ? 0 : length(hits(d))
Base.length(ray::AbstractRay; opl::Bool = false)
# Do
Base.length(d::Detector)                    # a detector is a hit container
optical_path_length(ray::AbstractRay)
isentering(r::AbstractRay) = isentering(r, intersection(r))
isentering(::AbstractRay, ::Nothing) = false
```
A verb whose meaning differs between types is not a rename: `set_pivot3d!(group, p)` moves only
the pivot, `position!(object, p)` moves the object.

**V2 · Small exported vocabulary, plain results.** Users get verbs, catalog constructors and
getters. The hooks a component author implements stay unexported, are documented on the
abstract type, and are `public` when dependent packages need them. Measurements return plain
tuples; a record read by field name is a small concrete struct (`PathSegment`), never a
`NamedTuple`. A mutating verb may return a useful value (`push!` returns the collection),
otherwise `nothing`.

**V3 · Constructors take intrinsic geometry and material only.** The object is created at the
origin along +y and placed afterwards with the kinematic verbs (beams and sources take position
and direction). One type parameter per physical argument, `T = promote_type(...)`; invariants are
checked here, once. Catalog variants are factory functions that return a generic type
(`RightAnglePrism(...) = Prism(RightAnglePrismSDF(...), n)`); a new struct only for a new
interaction or composition. Defining quantities are positional with SI defaults (`λ = 1e-6`),
options are keywords (`num_rings`, `M2`). Argument order: context, actor, carrier
(`interact3d(system, object, beam, ray)`).

### T · Types

**T1 · An abstract type is a contract.** Its docstring has `# Implementation reqs.` with
`## Fields`, `## Getters/setters`, `## Functions`. Default getters on the abstract type read
conventionally named fields (`pos`, `dir`, `shape`, `n`, `λ`, `beams`); siblings name the same
thing the same way, and code outside a type reads through the getter. Every interface verb has a
fallback on the abstract type, never a bare `function f end`: a missing required hook throws
`error(lazy"<verb> not defined for $(typeof(x))")`; `interact3d` without a method warns and
returns `nothing` (the object blocks the beam); a setter for irrelevant state returns `nothing`.

**T2 · Types carry information.** Fields are concrete: a type parameter for a field that would
otherwise be abstract, no `::Real`, `::Function`, `::Union{Nothing, Real}` or untyped `Vector`
fields, no parameters nobody dispatches on. Absence is `nothing` (`Nullable{T}`); a named "none"
singleton only as a value of an existing family that dispatch consumes (`NoBoundingSphere`,
`NoSampling`). Variants are singleton types, not `Symbol`s or `Val`s (exceptions: Makie-style
keywords such as `color = :wavelength`, and `hole_axis = :collimated`). Composites are objects
made of existing objects: `shape_trait_of(::C) = MultiShape()`, `interact3d` delegates to the
hit part and adds only the glue. Mutable for simulation state, immutable for wrappers; no
`_`-prefixed types.

### D · Dispatch and traits

**D1 · Every function is typed at least in its first argument.** Type further arguments as far as
unambiguous dispatch needs. Untyped parameters invite type piracy and ambiguities. Exceptions:
(1) the entry point of an opt-in trait API (D3); (2) a function as first argument (`do` syntax):
type the main object; (3) pure numeric formulas on numbers, vectors, points and matrices, so AD
and units pass (`rayleigh_range(λ, w0, M2, n = 1)`). An untyped argument is otherwise allowed
only if the method ignores it or passes it on unchecked.

**D2 · Choose behavior by dispatch.** No `isa`, `typeof` comparisons, `hasmethod`,
`haskey(kwargs, ...)` or Bool pseudo-traits to pick a code path. `isnothing` in plain control
flow and `isa` in input validation that throws are fine.

**D3 · Traits `<question>_of(x)` for questions across the hierarchy.** Trait values are empty
singletons without a `Trait` suffix; the public verb passes the trait as the first argument to
implementations typed on it. No trait for a question one verb asks, no predicate wrappers
(`is_static`). Within one hierarchy the default is a typed method on the root. An opt-in trait
across unrelated hierarchies (the kinematic API) throws for types that do not opt in, its entry
point is untyped and only forwards, and a type customizes at trait level. New verbs follow
[AbstractKinematicTrait.jl](src/AbstractTypes/AbstractKinematicTrait.jl): entry point → trait
method → derived methods built from public verbs → per-type primitives only where no verb
suffices (`rotate3d!` has four public signatures and one primitive, `rotate3d!(x, R)`).
```julia
kinematic_trait_of(::AbstractObject) = Movable(Oriented())
kinematic_trait_of(x) = error(lazy"kinematic_trait_of not defined for $(typeof(x))")
# A new verb point_at3d! (hypothetical)
point_at3d!(x, target) = point_at3d!(kinematic_trait_of(x), x, target)
point_at3d!(::Static, x, target) = static_error(x)
point_at3d!(::Movable, x, target::AbstractVector) = align3d!(x, target - position(x))
# Don't: kinematic_trait_of(x) = Static()   (silent default)
```

### S · State

**S1 · Small orchestration, explicit context, no hidden effects.** A driver loop does not branch
on component types and fits on one screen; complexity goes into small domain structs
(`Intersection`, hit records). Settings that change results come as arguments; global or
task-local state only as a cache that changes no result, behind existing verbs. `interact3d`
returns an interaction or `nothing`; only the solver mutates beams. A verb changes nothing beyond
its name, except data its own action invalidates (moving a detector empties its hits), and says
so. A precondition belongs to the function whose contract needs it, not to every caller.
Tolerances that change a physical result come from `Config`; algorithm constants stay inline with
a label.

### N · Naming

**N1 · Names follow roles.** Getter `name(x)`, setter `name!(x, new)`, adjacent one-liners; no
`get_`/`set_`/`calc_` for object state (except `Config` getters, backend helpers and
`set_pivot3d!`). Quantities are nouns (`electric_field`, `optical_power`); predicates `isx`
(`is_x_y` when clearer); `3d` only for operations in space (`translate3d!`, `intersect3d`). Short
field, long accessor (`pos` → `position`, `λ` → `wavelength`). Physics symbols for data (`λ`,
`θ`, `E0`, `w0`), English words for functions. Arguments by role (`object::AbstractObject`,
`ray`, `system`; concrete components by acronym, `cbs`), not `x`, `obj`, `sys`. Types
`Abstract<Role>`, catalog `[Shape][Variant][Family]`, backend suffix (`BoxSDF`). Constants
lowercase (`c_vacuum`), no SCREAMING_CASE; keyword families `num_*`, `*_max`, `show_*`.

### I · Internals

**I1 · The `_` prefix has three roles:** an internal protocol hook listed in the abstract type's
docstring, a private twin of a public verb (`_render!`), or a helper with three or more call sites
in the package. Any other unexported helper gets a plain name (`static_error`, `solve_beam!`); a
coherent multi-line step may be its own plain-named function so it can be tested; a one-line
helper with one call site stays inline.

**I2 · Factor shared logic, every time.** Logic shared by several models or parts becomes one
generic method, helper or loop over the parts, not parallel copies and not a `Symbol` table with
`getfield`. A check or formula that differs from an existing one only in a constant extends the
existing one (an argument or a method), it is not copied under a new name. One file per family: abstract type and contract, defaults, then each concrete type;
keyword-heavy signatures vertically.
```julia
# Do
for b in (beam.chief, beam.waist, beam.divergence)
    refractive_index!(rays(b)[id], n_new)
end
# Don't: a sibling check that differs only in its bound
check_grid_angle(θ) = 0 ≤ θ < π / 2 || throw(ArgumentError(...))   # next to check_half_angle(θ)
# Do: give the existing check the bound
check_half_angle(θ, θ_max = π)
```

### C · Comments and docs

**C1 · Comments label, docs explain.** A comment in a function body is a few capitalized words
without period (`# Catch first valid intersection`), dangers in capitals, gaps as `# FIXME`,
never `TODO`. Docstrings are short and correct: full ones on abstract types (T1) and user-facing
constructors (units such as `\\[m\\]` in the Julia source, orientation, limits as `!!! warning`,
`# Arguments`);
internal functions get one or two sentences. Known approximations are documented, not guarded on
hot paths. Explanations (why, flows) go on the docs page that embeds the docstrings: every core
system (a protocol or data flow across several types or files, e.g. the intersect-interact loop,
detector readout, sources, kinematics, the render protocol) has a Mermaid flowchart there in the
style of "Diagrams" in [docs/src/api/docdev.md](docs/src/api/docdev.md), derived from the code
and updated by every change to the system. Pages embed docstrings instead of repeating them; see
also "Documentation philosophy" and "Docstring conventions" in docdev.md.

### E, X · Errors and tests

**E1 · Exceptions by relevance:** `ArgumentError` for invalid input, `error(...)` for a missing
implementation, an own type only when callers catch it. Messages are short and name the type or
method; no `@assert` for input validation.

**X1 · Test physics through the public API** against analytic values or published designs with a
justified tolerance; move and rotate setups to show placement does not matter. Dummy subtypes
test the extension interface and fallbacks. Plain-named steps (I1) may be tested directly, `_`
functions not.

## Reviewing a pull request

A reviewing agent merges the branch with the current base, runs the test modules of the touched
code, and checks the diff against "Code style", reporting each violation with its rule ID,
location, and the expected behavior. A question the rules do not settle goes to the maintainers
as a question, not as a finding.

**Recording decisions:** when the maintainers settle a style question in a review, the PR that
applies the decision also adds it to "Code style", as a rule, an exception or an example, so no
agent has to ask again.

## Running Julia

- BMO requires Julia ≥ 1.12 (`Project.toml`). Use the newest Julia installed, not merely the
  minimum a Manifest allows. Machine-specific paths belong in an untracked file such as
  `CLAUDE.local.md`, not here.
- The repository is a Pkg workspace: `test` and `docs` are subprojects that pick up the local
  checkout automatically.

## Tests

- Single test module: `julia --project=test -e 'using Pkg; Pkg.instantiate()'` once, then
  `julia --project=test test/<path>.jl`. Run the modules of the code you change.
- Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'`, once before handing over.
- Test-only dependencies go into `test/Project.toml`, never the root `Project.toml`.
- New test files are `module TestXyz ... end`, included in [test/runtests.jl](test/runtests.jl).
  `Rendering/TestRenderErrors.jl` must run before anything loads the Makie extension.

## Agent skill

[skills/beamletoptics/](skills/beamletoptics/) teaches coding agents to *use* BMO. It ships with
every release, so it must describe exactly that release.

- **The author of a patch is responsible for keeping the skill correct.** A patch that changes
  public behavior (exports, constructor signatures, keyword arguments or defaults, units, sign or
  orientation conventions, validity limits of a beam model) updates the affected files
  (`SKILL.md`, `API.md`, `CONVENTIONS.md`, `CHECKLIST.md`, `components/`, `templates/`) in the
  same patch. Code in the skill must run as written: names that are `public` but not exported
  carry the prefix `BeamletOptics.`.
- [test/TestAgentSkill.jl](test/TestAgentSkill.jl) catches mechanical drift only (export table
  vs. `names(BeamletOptics)`, version field, `install_agent_skill`, headless templates). A changed
  meaning still passes; checking it is the author's job.
- On a minor version bump, review the whole skill and update the version field and the version
  named in the `SKILL.md` body.
- `templates/07_render_system.jl` needs GLMakie and is not in the test suite. After rendering
  changes, run `julia --project=docs skills/beamletoptics/scripts/run_templates.jl --render`.
- `BeamletOptics.install_agent_skill()` ([src/AgentSkill.jl](src/AgentSkill.jl)) copies
  `skills/beamletoptics/` out of the installed package; moving the folder means updating it and
  [docs/src/tutorials/ai_assistants.md](docs/src/tutorials/ai_assistants.md).

## Documentation build

[docs/src/api/docdev.md](docs/src/api/docdev.md) describes the build (DocumenterVitepress, output
in `docs/build/1`) and the figure pattern (scripts in `docs/src/assets`, loaded via
`Main.DocUtils.conditional_include`). Not obvious from that page:

- `GLOBAL_USE_PLACEHOLDERS` in [docs/DocUtils.jl](docs/DocUtils.jl) stays `true`; set it to
  `false` only to regenerate figures and back before committing. Such a run fills
  `docs/figure_cache` (untracked); delete a cached figure after changing its script.
- **GLMakie is the Makie backend** for figures, checks of `render!` and anything under `ext/`;
  on headless Linux run it under `xvfb-run -a`.
- **Windows link bug:** Documenter reads `[text]` followed by a parenthesized aside as a link;
  escape literal brackets as `\[1/m\]`.
- `@example` blocks end with `nothing # hide` or a trailing `;`, otherwise `show` output leaks
  into the page.
