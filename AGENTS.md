# BeamletOptics.jl (BMO): developer instructions for coding agents

A Julia package for non-sequential 3D ray and Gaussian beamlet tracing, used to simulate
breadboard optical setups (e.g. laser interferometers) with lenses, mirrors, beamsplitters,
polarizers and detectors. See [docs/src/index.md](docs/src/index.md) for the pitch and
[docs/src/api/conventions.md](docs/src/api/conventions.md) for the binding physical and
geometric conventions (global optical axis +y, right-handed frames, CCW rotations).

This file is for **developing** BMO. Guidance for **using** BMO lives in the agent skill
[skills/beamletoptics/](skills/beamletoptics/SKILL.md) (see "Agent skill" below).

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
`interact3d(system, object, beam, ray)` method (plus `intersect3d` if it needs custom
geometry), and the rest of the API (kinematics, threading, retracing) works without further
integration. A component may also add a `card_rows` method (and `card_actions`) to show its own rows
on its card in the interactive GUI, which lives in the separate package
[BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl) (recipe and developer
instructions there). When adding infrastructure,
prefer pushing complexity into the generic solver over asking component authors to handle it.

Current exception: coincident-boundary disambiguation (plate beamsplitters, cemented
doublets) is handled per component by returning a `Hint`
([src/AbstractTypes/AbstractSystem.jl](src/AbstractTypes/AbstractSystem.jl), consumed in
[src/System.jl](src/System.jl)). core.md explicitly calls this a burden on the developer.

**When reviewing or writing code, be suspicious of:**

- a component or algorithm that assumes a world-space orientation, a particular object order
  in a `System`, or that two objects are "adjacent" without deriving it from geometry at
  trace time;
- paraxial or near-normal-incidence shortcuts without an explicit, documented restriction on
  the angle of incidence;
- new functionality that makes users declare sequence or structure that the solver could
  infer from the scene.

## Code map

- `src/AbstractTypes/`: the interfaces (`AbstractObject`, `AbstractShape`, shape traits,
  `AbstractRay`/`Intersection`, `AbstractBeam`, `AbstractSystem`/`Hint`, kinematic trait).
- `src/System.jl`: the intersect-interact loop, tracing and retracing.
- `src/Rays.jl`, `PolarizedRays.jl`, `Beam.jl`, `Gaussian.jl`, `AstigmaticGaussian.jl`,
  `BeamGroups/`: beam models and sources.
- `src/OpticalComponents/`: components, one family per folder or file.
- `src/SDFs/`, `src/Mesh.jl`: geometry backends.
- `src/Exports.jl`: the public API. Changing it affects the agent skill (below).
- `ext/`: `BeamletOpticsMakieExt` and its `Render*.jl` files: `render!` and the live rendering
  (`RenderLive.jl` for objects and systems, the beam files for beams) behind the render handle
  protocol of `src/Render.jl`.
- The interactive GUI (`live_view`, cards, kinematic controls, view cube) is the separate package
  [BeamletOpticsGUI](https://github.com/StackEnjoyer/BeamletOpticsGUI.jl). It uses the exported
  names, the render handle protocol and the names declared `public` in `src/Exports.jl` (see
  "Public API for dependent packages" in [docs/src/api/api.md](docs/src/api/api.md)): renaming or
  changing any of these breaks the GUI, so treat them like exported API. It also uses some
  internal names (e.g. the abstract types, the kinematic and shape traits, `intersection`,
  `hits`, `objects`), which stay internal and may change: before renaming or changing an
  internal name, search the GUI for it and update the GUI in step.
- `skills/beamletoptics/`: the user-facing agent skill.

## Running Julia

- BMO requires Julia ≥ 1.12 (`Project.toml`). Use the newest Julia installed on the machine,
  not merely the minimum a Manifest allows. Machine-specific paths (e.g. the location of
  `julia.exe` when there is no `juliaup`) belong in an untracked local file such as
  `CLAUDE.local.md` (git-ignored), not here.
- The repository is a Pkg workspace: `test` and `docs` are subprojects that pick up the local
  checkout automatically.

## Tests

From [docs/src/api/contribute.md](docs/src/api/contribute.md):

- Single test module: `julia --project=test -e 'using Pkg; Pkg.instantiate()'` once, then
  `julia --project=test test/<path>.jl`.
- Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'` (or `] test`).
- Test-only dependencies go into `test/Project.toml`, never the root `Project.toml`.
- New test files are `module TestXyz ... end` and must be included in
  [test/runtests.jl](test/runtests.jl). Order matters: `Rendering/TestRenderErrors.jl` must run
  before anything loads the Makie extension.

## Agent skill

[skills/beamletoptics/](skills/beamletoptics/) teaches coding agents to *use* BMO. It ships
with every release, so it must describe exactly that release.

- **The author of a patch (usually a coding agent) is responsible for keeping the skill correct
  for the release the patch lands in.** A patch that changes public behavior (exports,
  constructor signatures, keyword arguments or their defaults, units, sign or orientation
  conventions, the validity limits of a beam model) reviews the affected skill files
  (`SKILL.md`, `API.md`, `CONVENTIONS.md`, `CHECKLIST.md`, `components/`, `templates/`) and
  updates them in the same patch. State in the PR description which skill files were updated,
  or that none were affected.
- [test/TestAgentSkill.jl](test/TestAgentSkill.jl) catches mechanical drift only. It fails if
  the export table in `API.md` differs from `names(BeamletOptics)`, if
  `metadata.beamletoptics-version` in `SKILL.md` differs from the major.minor version in
  `Project.toml`, if `install_agent_skill` breaks, or if a headless template (01–06) throws.
  Passing tests do not show that the skill is correct: a changed meaning, e.g. a flipped sign
  convention, still runs. Checking that is the author's job (above).
- On a minor version bump, review the whole skill and update both the version field and the
  version named in the `SKILL.md` body.
- `templates/07_render_system.jl` needs GLMakie and is not part of the test suite. After
  rendering changes, run all templates with
  `julia --project=docs skills/beamletoptics/scripts/run_templates.jl --render`.
- Users install the version-matched copy with `BeamletOptics.install_agent_skill()`
  ([src/AgentSkill.jl](src/AgentSkill.jl)), which copies `skills/beamletoptics/` out of the
  installed package. Keep the folder at that path, or update the function and
  [docs/src/tutorials/ai_assistants.md](docs/src/tutorials/ai_assistants.md) with it.

## Docstrings

Follow "Documentation philosophy" and "Docstring conventions" in
[docs/src/api/docdev.md](docs/src/api/docdev.md). In short: pages embed docstrings instead
of repeating them, constructor docstrings must be self-sufficient in the REPL (units,
conventions, limitations), and abstract type docstrings state the interface a subtype
implements.

## Documentation build

[docs/src/api/docdev.md](docs/src/api/docdev.md) describes the build (DocumenterVitepress,
output in `docs/build/1`, served with LiveServer) and the figure pattern (scripts in
`docs/src/assets`, loaded via `Main.DocUtils.conditional_include`). Not obvious from that page:

- `GLOBAL_USE_PLACEHOLDERS` at the top of [docs/DocUtils.jl](docs/DocUtils.jl) switches
  local builds between real figures and fast placeholders. Keep it `true`; set it to `false`
  only to regenerate figures, and set it back before committing. CI always renders. Such a run
  fills `docs/figure_cache` (untracked), which later placeholder builds use instead of
  placeholders; delete a cached figure after changing its script.
- **GLMakie is the Makie backend.** Use it for figures, ad-hoc checks of `render!` and anything
  under `ext/`. On headless Linux, run it under `xvfb-run -a` (CI does the same).
- **Windows link bug:** Documenter reads `[text]` followed by a parenthesized aside, e.g.
  `` `α` in [1/m] (Lambert-Beer: ...) ``, as a link and fails with "colons not allowed in
  paths". Escape literal brackets as `\[1/m\]`.
- `@example` blocks must end with a suppressing statement (`nothing # hide` or a trailing
  `;`), otherwise `Base.show` output of the last expression leaks into the page.

## Code conventions

- Style baseline: the [SciML Style Guide](https://github.com/SciML/SciMLStyle) (not strictly
  enforced).
- **Preconditions belong to the function whose contract needs them**, not to every call site
  that reaches it. If a public function must establish some state (e.g. a fresh trace), it
  does so itself, and callers one layer up trust it instead of repeating the setup.
- **Extend verbs by dispatch, not by `_helper` functions.** Before writing an unexported
  `_helper`, check whether the subtask is one of the following, each of which has a dispatch
  answer:
  - *Normalizing an argument* (point or object, axis and angle or matrix, with or without
    pivot): add a method of the same verb that converts the argument and calls the verb
    again, e.g. `rotate3d!(x, axis, θ)` → `rotate3d!(x, R)`.
  - *Branching on a type or trait* (`isa`, `isnothing(kw)` selecting a code path): dispatch
    on the type or on a trait (`kinematic_trait_of`, `shape_trait_of`). An optional argument
    that changes behavior is a positional method, not a `kw = nothing` keyword.
  - *Repeating an existing operation* (rotating about a pivot, aligning a direction): call
    the existing public verb. Derived methods call public entry points only, so type-specific
    methods of those verbs apply automatically.
  - *Checking a precondition a trait already encodes* (static, directed vs. oriented): let
    the trait branch handle it.

  New verbs follow the pattern of [AbstractKinematicTrait.jl](src/AbstractTypes/AbstractKinematicTrait.jl):
  public entry point → trait method → derived methods built from existing verbs →
  per-type primitives only where no existing verb suffices. `rotate3d!` has four public
  signatures but only one per-type primitive, `rotate3d!(x, R)`:

  ```julia
  point_at3d!(x, target) = point_at3d!(kinematic_trait_of(x), x, target)
  point_at3d!(::Static, x, target) = _static_error(x)
  point_at3d!(t::Movable, x, target) = point_at3d!(t, x, position(target))   # any positioned target
  point_at3d!(::Movable, x, target::AbstractVector) = align3d!(x, target - position(x))
  ```

  A `_helper` is fine for type-independent work that has no natural verb: error
  constructors (`_static_error`), numeric kernels, internal data plumbing. It must not
  branch on argument types.
- New or changed functionality comes with tests, docstrings, and docs or examples where
  relevant ([docs/src/api/contribute.md](docs/src/api/contribute.md)).
