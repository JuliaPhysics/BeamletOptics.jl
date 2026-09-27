# BeamletOptics Agent Skill

This folder is an [Agent Skill](https://docs.claude.com/en/docs/agents-and-tools/agent-skills/overview)
that teaches coding agents how to write simulations with
[BeamletOptics.jl](https://github.com/JuliaPhysics/BeamletOptics.jl).

Install it in one of these ways:

- **Recommended:** copy the version that matches your installed package from Julia:
  `using BeamletOptics; BeamletOptics.install_agent_skill()` installs to `./.claude/skills/beamletoptics`.
  Pass `joinpath(homedir(), ".claude", "skills")` for a personal installation. Rerun after `Pkg.update`.
- `npx skills add JuliaPhysics/BeamletOptics.jl` (the CLI discovers `skills/beamletoptics/SKILL.md`).
  This copy follows the development branch and may be newer than your installed release.
- Copy this folder manually to `~/.claude/skills/beamletoptics/` (personal) or
  `<project>/.claude/skills/beamletoptics/` (project).

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
BeamletOptics installed. `07_render_system.jl` additionally needs `GLMakie` and a display
(`xvfb-run -a` on headless Linux).

Paths below are relative to this skill directory; `<env>` is a Julia project that has BeamletOptics
(and `GLMakie` for `--render`) installed:

```sh
julia --project=<env> templates/03_psf_airy.jl
julia --project=<env> scripts/run_templates.jl            # smoke-test 01-06
julia --project=<env> scripts/run_templates.jl --render   # all templates, needs GLMakie
```

From the root of the BeamletOptics repository, use the package env (templates 01 to 06) or the docs
env (all templates):

```sh
julia --project=. skills/beamletoptics/scripts/run_templates.jl
julia --project=docs skills/beamletoptics/scripts/run_templates.jl --render
```

## Helper scripts

- `scripts/api_lookup.jl NAME...`: print docstrings and method signatures (`--exports`, `--subtypes TYPE`)
- `scripts/run_templates.jl`: run every template and report failures (useful after API changes)

## Maintenance

When the public API changes (`src/Exports.jl`, constructor signatures, keyword arguments,
conventions), update the matching file in `components/` and rerun the templates (see above).
Developer instructions for the package itself are in `AGENTS.md` at the repository root.

`test/TestAgentSkill.jl` guards two things in the package test suite:

- the exported-names table in `API.md` must match `names(BeamletOptics)` in both directions;
- `metadata.beamletoptics-version` in `SKILL.md` must match the major.minor version in `Project.toml`.
  When bumping the minor version, review the skill and update this field (and the version named in
  the `SKILL.md` body).
