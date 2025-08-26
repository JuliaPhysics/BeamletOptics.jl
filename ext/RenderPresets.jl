render!(ax::_RenderEnv, refr::BMO.AbstractRefractiveOptic; kwargs...) = _render!(ax, refr; transparency=true, color=:white, kwargs...)
render!(ax::_RenderEnv, refl::BMO.AbstractReflectiveOptic; kwargs...) = _render!(ax, refl; transparency=false, color=:silver, kwargs...)

render!(ax::_RenderEnv, lens::Lens; kwargs...) = _render!(ax, lens; transparency=true, color=RGBf(0.678, 0.847, 0.902), alpha=0.5, kwargs...)

render!(ax::_RenderEnv, bs::ThinBeamsplitter; kwargs...) = _render!(ax, bs; transparency=true, color=:magenta, kwargs...)

render!(ax::_RenderEnv, nino::BMO.NonInteractableObject; kwargs...) = _render!(ax, nino; transparency=false, color=:grey, kwargs...)

render!(ax::_RenderEnv, nino::BMO.IntersectableObject; kwargs...) = _render!(ax, nino; transparency=true, color=:grey, alpha=0.1)