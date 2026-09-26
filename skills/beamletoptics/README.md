# BeamletOptics Agent Skill

This folder is an [Agent Skill](https://docs.claude.com/en/docs/agents-and-tools/agent-skills/overview)
that teaches coding agents how to write simulations with
[BeamletOptics.jl](https://github.com/JuliaPhysics/BeamletOptics.jl).

Install it as:

- `npx skills add JuliaPhysics/BeamletOptics.jl` (the CLI discovers `skills/beamletoptics/SKILL.md`), or
- copy this folder to `~/.claude/skills/beamletoptics/` (personal) or `<project>/.claude/skills/beamletoptics/` (project).

The canonical entry point is `SKILL.md`.

## Contents

- `SKILL.md`: main skill definition (frontmatter and instructions)
- `API.md`: compact API primer (build, solve, evaluate)
- `CONVENTIONS.md`: units, coordinate system, sign conventions
- `WORKFLOW.md`: recommended patterns (scans, focus search, interferometers, groups)
- `VISUALIZATION.md`: rendering with the Makie extension
- `CHECKLIST.md`: pre-delivery checklist and common pitfalls
- `components/`: per-component reference (constructors, orientation, caveats)
- `templates/`: runnable example scripts
- `scripts/`: helper scripts

## Templates

The files in `templates/` are small, complete scripts. Templates 01 to 06 run headless with only
BeamletOptics installed. `07_render_system.jl` additionally needs `CairoMakie`.

```sh
julia --project=<env with BeamletOptics> templates/03_psf_airy.jl
julia --project=<env with BeamletOptics> scripts/run_templates.jl            # smoke-test 01-06
julia --project=<env with BeamletOptics + CairoMakie> scripts/run_templates.jl --render
```

Inside this repository, `--project=.` (package env) runs 01 to 06 and `--project=docs` runs all of them.

## Helper scripts

- `scripts/api_lookup.jl NAME...`: print docstrings and method signatures (`--exports`, `--subtypes TYPE`)
- `scripts/run_templates.jl`: run every template and report failures (useful after API changes)

## Maintenance

When the public API changes (`src/Exports.jl`, constructor signatures), update the matching file in
`components/` and rerun `scripts/run_templates.jl --render`.
