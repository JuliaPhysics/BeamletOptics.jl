render!(ax::_RenderEnv, refr::BMO.AbstractRefractiveOptic; kwargs...) = _render!(ax, refr; transparency=true, color=:white, kwargs...)
render!(ax::_RenderEnv, refl::BMO.AbstractReflectiveOptic; kwargs...) = _render!(ax, refl; transparency=false, color=:silver, kwargs...)

render!(ax::_RenderEnv, bs::ThinBeamsplitter; kwargs...) = _render!(ax, bs; transparency=false, color=:magenta, kwargs...)

render!(ax::_RenderEnv, nino::BMO.NonInteractableObject; kwargs...) = _render!(ax, nino; transparency=false, color=:grey, kwargs...)