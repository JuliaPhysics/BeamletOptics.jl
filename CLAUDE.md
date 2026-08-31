# BeamletOptics.jl (BMO)

A Julia package for 3D Gaussian beamlet tracing, used to simulate breadboard optical
setups (e.g. laser interferometers) with components like lenses, mirrors, beamsplitters
and coatings. See [docs/src/index.md](docs/src/index.md) for the full pitch and
[docs/src/api/conventions.md](docs/src/api/conventions.md) for binding physical/geometric
conventions (global optical axis = +y, right-handed coordinate systems, CCW rotations).

## Design philosophy

BMO is meant to feel like a **digital laboratory**: the user places components in
3D space the way they would physically arrange them on an optical breadboard, and the
tracer figures out everything else on its own. This is the lens through which every API
and architecture decision should be evaluated. Concretely, from
[docs/src/api/core.md](docs/src/api/core.md):

1. Optical interaction is decoupled from geometry representation.
2. Optical elements are closed volumes (or must mimic being one) — exceptions apply,
   e.g. coatings.
3. Elements must be freely movable and work correctly for (almost) any angle of
   incidence — no paraxial shortcuts, no assumed canonical orientation.
4. Without extra knowledge, tracing is non-sequential: the solver discovers what a
   ray/beam hits next by searching the scene, not by a user-declared optical path or
   object order.
5. With extra knowledge (a `Hint`), tracing can go sequential for performance/ambiguity
   resolution — but this is an optimization, never a requirement placed on the user.

**The extension promise:** a developer defines a new `AbstractObject` subtype and its
`interact3d(system, object, beam, ray) -> AbstractInteraction` method (plus `intersect3d`
if it needs custom geometry) — and the rest of the API (kinematics, threading, the `Hint`
mechanism, ...) works without further integration work. When adding infrastructure, prefer
pushing complexity into the generic solver over asking component authors to handle it (see
the `Hint` interface in [System.jl](src/System.jl) as the reference example of that
trade-off). Note: coincident-boundary disambiguation (coatings, cemented doublets) is
currently still handled per-component via `Hint` + shape-identity comparison, not
centrally by the solver — `MultiIntersection` (`src/AbstractTypes/AbstractIntersection.jl`)
is scaffolding towards changing that.

**Corollary — when reviewing or writing code, be suspicious of:**
- Any component or algorithm that assumes a specific world-space orientation, a
  particular object ordering in a `System`, or that two objects are "adjacent" without
  deriving that from geometry at trace time.
- Paraxial or near-normal-incidence shortcuts introduced without an explicit,
  documented restriction on angle of incidence.
- New functionality that requires users to manually declare sequence/structure that the
  solver could in principle infer from the scene.

## Running Julia

Always use the newest Julia version installed on the machine to run or build BMO, not
just whatever a pinned Manifest specifies as a minimum. There is no `juliaup` on this
machine, so pick the right `julia.exe` explicitly, e.g.:

```
"/c/Users/uitt_hu/AppData/Local/Programs/Julia-1.12.1/bin/julia.exe" --project=. -e '...'
```

Re-check installed versions each session (`AppData\Local\Programs\Julia-<version>\bin`)
rather than assuming a specific version stays the newest.

## Documentation philosophy & docstring style

The full, canonical guide now lives in
[docs/src/api/docdev.md](docs/src/api/docdev.md) ("Documentation philosophy" and
"Docstring conventions" sections) — keep that page and this summary in sync when either
changes. Key points:

- **Docs live off docstrings and figures, not hand-duplicated prose.** A `.md` page
  should embed docstrings via `@docs` blocks and Makie-rendered figures rather than
  re-explain in free text what a docstring/figure should already cover. If a page needs
  to explain something a docstring is missing, extend the docstring — don't compensate
  for it in the page.
- **Standard structure:** indented signature-synopsis line, blank line, prose
  description using `[`Type`](@ref)`/`[`func`](@ref)` cross-references, then an optional
  `# Fields` (structs) or `# Arguments` (functions) bullet list (`- `name`: description`).
  Extra `# <Section>` headers (e.g. `# Sign convention`) are fine if an existing one
  doesn't already fit; don't repeat prose and bullets.
- **Constructor docstrings must be self-sufficient for the REPL.** `?MyComponent` is
  often a user's only source of information — the constructor docstring must cover units,
  sign/orientation conventions, physical assumptions and limitations, not just parameter
  types.
- **Type/abstract type docstrings must state the interface they define or fulfill** —
  which functions a subtype needs to implement, or which trait/interface a concrete type
  participates in. Reference pattern: `SingleShape`/`MultiShape` in
  `src/AbstractTypes/AbstractShapeTrait.jl` (an `# AbstractObject implementation reqs.`
  section).

## Code conventions

- **DRY: precondition/invariant-establishing calls live in the function whose contract
  needs them, not at every call site that happens to reach it.** If a function's contract
  requires some setup (e.g. `_reset_beam!` establishing "always trace fresh" for
  `trace_system!`), that call belongs inside the contract-owning function itself — not
  duplicated by callers that merely dispatch into it. A caller one layer up should trust
  the callee to uphold its own contract rather than re-establishing the precondition
  itself; doing so is dead work, not useful defensive redundancy. Reference case:
  `trace_system!` (all methods, `src/System.jl`) calls `_reset_beam!` as its first
  statement because it is public API also called directly (tests, docs), so it must
  guarantee a fresh trace regardless of caller. `solve_system!` used to also call
  `_reset_beam!` right before dispatching into `trace_system!` via `solve_leaf!` — that
  outer call was redundant (nothing ran in between to justify re-establishing the same
  state) and was removed rather than centralized upward, since moving it to
  `solve_system!` only and stripping it from `trace_system!` would have broken
  `trace_system!`'s standalone contract (see `test/TestSystem.jl`, which calls
  `trace_system!` directly twice on the same beam and expects each call to retrace
  fresh).

## Documentation build

Building the docs locally is described in
[docs/src/api/docdev.md](docs/src/api/docdev.md) — `] dev BeamletOptics` inside the `docs`
environment, then run `docs/make.jl`. Two things that aren't obvious from that page:

- **`docs/DocUtils.jl`**: `GLOBAL_USE_PLACEHOLDERS` (top of file) controls whether local
  builds render real figures via GLMakie/CairoMakie or fast placeholder images. Keep it
  `true` for day-to-day local builds — a full render pass takes several minutes. Only
  flip it to `false` temporarily when you actually need to regenerate real figures, then
  set it back before committing.
- **GLMakie is BMO's assumed/default Makie backend.** CairoMakie is only used where a
  static `.png` export is specifically wanted (e.g. some doc figures). When manually
  exercising `render!` or anything under `ext/` (`BeamletOpticsMakieExt` and its
  constituent `Render*.jl` files) — e.g. an ad-hoc smoke test after touching rendering
  code — load GLMakie, not CairoMakie, unless a static image is explicitly what's needed.
- **Windows-only link bug**: Documenter treats a bare `[text]` followed by a
  parenthesized aside — e.g. `` `α` in [1/m] (Lambert-Beer: ...) `` — as a malformed
  markdown link. On Windows this fails the build with "colons not allowed in paths"
  because the aside text becomes the bogus link target. Escape literal square brackets in
  docstrings/markdown as `\[1/m\]` instead of `[1/m]` when they aren't meant to be a link
  or footnote reference.
- Example-generating scripts in `@example` blocks must end with a suppressing statement
  (`nothing # hide` or a trailing `;`) — otherwise the last expression's `Base.show`
  output (e.g. a full `System`/`Lens` struct dump) leaks into the rendered page and the
  build console.

## Contributing expectations

Per [docs/src/api/contribute.md](docs/src/api/contribute.md): new/changed functionality
should come with tests, docstrings, and documentation/examples where relevant. Style
baseline is the [SciML Style Guide](https://github.com/SciML/SciMLStyle) (not strictly
enforced).
