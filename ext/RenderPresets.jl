# Coating visualization palette
coating_color(::BMO.Uncoated, default_color = RGBAf(0.678, 0.847, 0.902, 0.4)) = default_color
coating_color(::BMO.SimpleARCoating, default_color = nothing) = RGBAf(0.18, 0.80, 0.44, 0.65)
coating_color(::BMO.SimpleHRCoating, default_color = nothing) = RGBAf(0.95, 0.77, 0.20, 0.95)
coating_color(::BMO.SimpleBeamsplitterCoating, default_color = nothing) = RGBAf(0.90, 0.49, 0.13, 0.70)
coating_color(::BMO.JonesCoating, default_color = nothing) = RGBAf(0.60, 0.20, 0.80, 0.70)
coating_color(::BMO.ThinFilmCoating, default_color = nothing) = RGBAf(0.10, 0.74, 0.61, 0.65)
coating_color(::BMO.AbstractCoatingModel, default_color = RGBAf(0.85, 0.45, 0.85, 0.70)) = default_color

render!(ax::_RenderEnv, refr::BMO.AbstractRefractiveOptic; kwargs...) = _render!(ax, refr; transparency=true, color=RGBf(0.678, 0.847, 0.902), alpha=0.5, kwargs...)
render!(ax::_RenderEnv, refl::BMO.AbstractReflectiveOptic; kwargs...) = _render!(ax, refl; transparency=false, color=:silver, kwargs...)

render!(ax::_RenderEnv, lens::Lens; kwargs...) = _render!(ax, lens; transparency=true, color=RGBf(0.678, 0.847, 0.902), alpha=0.5, kwargs...)
render!(ax::_RenderEnv, lens::DoubletLens; kwargs...) = _render!(ax, lens; transparency=true, color=RGBf(0.678, 0.847, 0.902), alpha=0.5, kwargs...)

render!(ax::_RenderEnv, c::Coating; kwargs...) = _render!(ax, c; transparency=true, color=coating_color(c.model), kwargs...)

render!(ax::_RenderEnv, nino::BMO.NonInteractableObject; kwargs...) = _render!(ax, nino; transparency=false, color=:grey, kwargs...)

render!(ax::_RenderEnv, nino::BMO.IntersectableObject; kwargs...) = _render!(ax, nino; transparency=true, color=:grey, kwargs...)