#=
Coloring of beam plots: in the color of the wavelength of each ray, `color = :wavelength` or
`color = (:wavelength, alpha)`, or in a single color. A live handle draws both with one color per
vertex, such that `render_settings!` switches between them without creating new plots.
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
    _SingleColor

The color of a beam plot that is the same for all of its vertices. Created from the `color` kwarg
by [`_coloring`](@ref).
"""
struct _SingleColor
    color::RGBAf
end

"""
    _coloring(color)

Interprets the `color` kwarg of a beam plot: `:wavelength` and `(:wavelength, alpha)` become a
[`_ByWavelength`](@ref), a single color that `Makie.to_color` knows, e.g. `:blue` or `(:blue, 0.3)`,
a [`_SingleColor`](@ref). Any other `color`, e.g. a vector of colors, is returned as it is and
passed on to `Makie`.
"""
_coloring(color) = _single_color(color)
_coloring(color::Symbol) = color === :wavelength ? _ByWavelength(1.0f0) : _single_color(color)
function _coloring(color::Tuple{Symbol, Real})
    first(color) === :wavelength || return _single_color(color)
    return _ByWavelength(Float32(last(color)))
end

function _single_color(color)
    c = try
        Makie.to_color(color)
    catch e
        e isa InterruptException && rethrow()
        # no color: `Makie` reports it when the plot is created
        return color
    end
    return c isa Makie.Colorant ? _SingleColor(RGBAf(c)) : color
end

"""
    _coloring(color, live::Bool)

The coloring with which a beam plot draws the `color`, see [`_coloring`](@ref). Of a plot that is
not `live`, i.e. one of `render!`, a single color is passed on to `Makie` as it is.
"""
function _coloring(color, live::Bool)
    c = _coloring(color)
    return (live || c isa _ByWavelength) ? c : color
end

"""
    _vertex_colors(coloring, n, by_wavelength) -> Union{Vector{RGBAf}, Nothing}

The colors of the `n` vertices of a beam plot drawn with the `coloring`: `by_wavelength(coloring)`
for a [`_ByWavelength`](@ref), `n` times the color of a [`_SingleColor`](@ref), and `nothing` for
any other `coloring`, which the plot keeps as its `color` attribute.
"""
_vertex_colors(c::_ByWavelength, n, by_wavelength) = by_wavelength(c)
_vertex_colors(c::_SingleColor, n, by_wavelength) = fill(c.color, n)
_vertex_colors(_, n, by_wavelength) = nothing

"""The display color of the wavelength `λ` [m] with the opacity `alpha`."""
_wavelength_rgba(λ, alpha::Real) = RGBAf(BMO.wavelength_color(λ)..., alpha)
_wavelength_rgba(λ, c::_ByWavelength) = _wavelength_rgba(λ, c.alpha)
