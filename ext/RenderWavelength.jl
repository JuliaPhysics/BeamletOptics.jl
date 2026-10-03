#=
Drawing beams in the color of their wavelength, `color = :wavelength` or `color = (:wavelength, alpha)`
=#

"""
    _ByWavelength

The color of a beam plot that is taken from the wavelength of each ray (or beamlet), see
[`wavelength_color`](@ref), with the opacity `alpha`. Created from the `color` kwarg by
[`_coloring`](@ref).
"""
struct _ByWavelength
    alpha::Float32
end

"""
    _coloring(color)

Interprets the `color` kwarg of a beam plot: `:wavelength` and `(:wavelength, alpha)` become a
[`_ByWavelength`](@ref), any other `color` is returned as it is and passed on to `Makie`.
"""
_coloring(color) = color
_coloring(color::Symbol) = color === :wavelength ? _ByWavelength(1.0f0) : color
function _coloring(color::Tuple{Symbol, Real})
    first(color) === :wavelength || return color
    return _ByWavelength(Float32(last(color)))
end

"""The display color of the wavelength `λ` [m] with the opacity `alpha`."""
_wavelength_rgba(λ, alpha::Real) = RGBAf(BMO.wavelength_color(λ)..., alpha)
_wavelength_rgba(λ, c::_ByWavelength) = _wavelength_rgba(λ, c.alpha)
