"""
    PolarizingIsolator

Simple optical isolator built from a `PolarizingCubeBeamsplitter` and a
`HalfWaveplate`. The waveplate is placed in front of the cube beamsplitter
and can be rotated to tune the polarization. Incoming horizontally polarized
light is transmitted while light traveling in reverse is reflected.
"""
struct PolarizingIsolator{T,B<:PolarizingCubeBeamsplitter{T},W<:Waveplate{T}}
    beamsplitter::B
    waveplate::W
end

shape_trait_of(::PolarizingIsolator) = MultiShape()
shape(iso::PolarizingIsolator) = (iso.beamsplitter, iso.waveplate)

"""
    PolarizingIsolator(size, n; offset=0)

Construct a polarizing isolator using a cube beamsplitter of side length `size`
and refractive index function `n`. The half-wave plate is positioned `offset`
meters in front of the cube.
"""
function PolarizingIsolator(size::Real, n::RefractiveIndex; offset::Real=0)
    bs = PolarizingCubeBeamsplitter(size, n)
    wp = HalfWaveplate(size)
    translate3d!(wp, [0, offset, 0])
    return PolarizingIsolator(bs, wp)
end

