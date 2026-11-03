# Using BMO with AI assistants

BMO ships an [agent skill](https://docs.claude.com/en/docs/agents-and-tools/agent-skills/overview): a folder of instructions, API references and runnable example scripts that teaches AI coding assistants such as Claude Code how to write BMO simulations. It covers units and coordinate conventions, the choice of beam model, all component constructors, detector readout, rendering and common pitfalls.

## Installation

The skill is part of the package, in the `skills/beamletoptics` folder of the installed package directory. Copy the version that matches your installed BMO into your project with

```julia
using BeamletOptics
BeamletOptics.install_agent_skill()   # -> ./.claude/skills/beamletoptics
```

For a personal installation that applies to all projects, pass the user-level skill directory:

```julia
BeamletOptics.install_agent_skill(joinpath(homedir(), ".claude", "skills"))
```

Run the command again after `Pkg.update` to keep the skill in sync with the package version.

```@docs; canonical=false
BeamletOptics.install_agent_skill
```

!!! info "Installing from GitHub"
    The skill can also be installed with `npx skills add JuliaPhysics/BeamletOptics.jl`. This copy follows the development branch and can be newer than your installed BMO release. The skill records the BMO version it describes (`beamletoptics-version` in `SKILL.md`) and instructs the assistant to verify signatures when the versions differ.

!!! tip "Developing BMO with an assistant"
    The skill describes how to *use* BMO. To change the package itself with an assistant, see [Developing BMO with an assistant](@ref) on the [Contributing](@ref) page.
