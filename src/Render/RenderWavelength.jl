# Color of a wavelength, used by ext/RenderWavelength.jl

"""
    wavelength_color(λ) -> NTuple{3, Float64}

The display color `(r, g, b)` (each in `[0, 1]`, sRGB) of light of the vacuum wavelength `λ` [m], e.g.
to draw rays in the color of their wavelength, see the `color = :wavelength` option of
[`render!`](@ref). Uses the piecewise linear approximation of Dan Bruton (violet at 380 nm over
blue, cyan, green, yellow and orange to red at 780 nm) with a gamma of 0.8, and the intensity falls
off to 30 % towards both ends of the visible range (380 to 420 nm and 700 to 780 nm).

Outside of 380 to 780 nm there is no visible color: the color of the nearest end of the spectrum
is returned, i.e. a dim violet for UV and a dim red for IR, such that rays of an invisible
wavelength are still drawn, and the color stays continuous. The result is a perceptual
approximation for plots, not a colorimetric conversion.

The function is plain Julia and needs no plotting backend. Throws an `ArgumentError` unless `λ` is
positive and finite.

```julia
wavelength_color(450e-9)   # blue
wavelength_color(650e-9)   # red
```
"""
function wavelength_color(λ::Real)
    (isfinite(λ) && λ > 0) || throw(ArgumentError("the wavelength must be positive and finite, got $λ"))
    # [nm], clamped to the visible range
    nm = clamp(1e9 * Float64(λ), 380.0, 780.0)
    r, g, b = if nm < 440
        (-(nm - 440) / 60, 0.0, 1.0)
    elseif nm < 490
        (0.0, (nm - 440) / 50, 1.0)
    elseif nm < 510
        (0.0, 1.0, -(nm - 510) / 20)
    elseif nm < 580
        ((nm - 510) / 70, 1.0, 0.0)
    elseif nm < 645
        (1.0, -(nm - 645) / 65, 0.0)
    else
        (1.0, 0.0, 0.0)
    end
    # intensity fall-off at the ends of the spectrum
    f = if nm < 420
        0.3 + 0.7 * (nm - 380) / 40
    elseif nm <= 700
        1.0
    else
        0.3 + 0.7 * (780 - nm) / 80
    end
    γ = 0.8
    return map(c -> clamp((f * c)^γ, 0.0, 1.0), (r, g, b))
end
