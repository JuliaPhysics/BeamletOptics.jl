"""
    OpenRaman

Custom `BeamletOptics` components for the [OpenRAMAN spectroscopy](@ref) tutorial.

This module is **not** a registered package. It is bundled with the BeamletOptics
documentation as a worked example of how to extend the package with your own optical
elements. It is loaded once from `docs/make.jl` and consumed both by
`openraman_showcase.jl` and by the tutorial page itself.

Two components are provided:

- [`ReflectiveGrating`](@ref): a multi-shape reflective diffraction grating implementing
  the vector grating equation
- [`LongpassDichroicMirror`](@ref): a dichroic mirror whose interaction branches on the
  wavelength of the incoming ray

alongside a small set of `SellmeierEquation` glass models used by the OpenRAMAN optics.
"""
module OpenRaman

using BeamletOptics
using LinearAlgebra: dot, normalize
# LScene is only needed to type the `render!` signatures. Without it, `render!(ax, ...)`
# would be ambiguous against the `render!(ax::_RenderEnv, ::AbstractObject)` fallback of
# the BeamletOptics Makie extension.
using GLMakie: LScene

const BMO = BeamletOptics

include("glasses.jl")
include("grating.jl")
include("longpass.jl")

export ReflectiveGrating, RectangularReflectiveGrating, LongpassDichroicMirror
export N_BK7, N_SF10, N_SF6HT, N_BAF10, UV_Fused_Silica

end # module OpenRaman
