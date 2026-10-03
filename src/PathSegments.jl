"""
    PathSegment{T}

Element type of the vector returned by [`path_segments`](@ref), a `NamedTuple` with the fields

- `start`, `stop`: start and end point of the segment, `Point3{T}` in [m]
- `s_start`, `s_stop`: geometric path length from the source to the start and the end of the segment [m]
- `opl_start`, `opl_stop`: optical path length from the source to the start and the end of the segment [m]
- `λ`: wavelength of the segment in [m]
- `depth`: depth of the beam in the beam tree, `1` for a root beam, `2` for its children, ...
- `parent`: index of the preceding segment in the vector (the previous segment of the same beam or,
  for the first segment of a child beam, the last segment of its parent beam), `0` for the first
  segment of a root beam
- `beam`: index of the root beam in a beam group, `1` for a single beam
- `final`: `true` if the segment is a final ray without intersection, drawn with the length `flen`
"""
const PathSegment{T} = @NamedTuple{
    start::Point3{T}, stop::Point3{T},
    s_start::T, s_stop::T, opl_start::T, opl_stop::T,
    λ::T, depth::Int, parent::Int, beam::Int, final::Bool
}

"""
    _chief_beam(beam::AbstractBeam)

The [`Beam`](@ref) of `beam` that carries its geometric path: the beam itself, or the chief beam of a
beamlet.
"""
_chief_beam(beam::Beam) = beam
_chief_beam(gauss::GaussianBeamlet) = gauss.chief
_chief_beam(agb::AstigmaticGaussianBeamlet) = agb.c

"""
    path_segments(beam::AbstractBeam; flen = 1.0) -> Vector{PathSegment}
    path_segments(group::AbstractBeamGroup; flen = 1.0) -> Vector{PathSegment}

Flattens the traced `beam` and all of its child beams (e.g. after beam splitting) into a vector of
straight segments, one per ray, in depth-first order. For a [`GaussianBeamlet`](@ref) or an
[`AstigmaticGaussianBeamlet`](@ref) the segments are those of the chief ray. For a beam group the
segments of all its beams are concatenated, the `beam` field gives the index of the beam.

Each segment is a [`BeamletOptics.PathSegment`](@ref) `NamedTuple` with the start and end point, the
accumulated geometric path length `s_start`, `s_stop` and optical path length `opl_start`,
`opl_stop` (refractive index of the medium of the ray times its length), the wavelength `λ`, the
`depth` in the beam tree and the `parent` segment index, from which the branches can be rebuilt. A
child beam continues from the end of its parent beam, i.e. its first segment starts at the path
length at which its parent ends. A final ray without an intersection has no length of its own: it
is `flen` [m] long and flagged by `final = true`.

!!! note
    The beam must have been traced with [`solve_system!`](@ref); an untraced beam consists of a
    single final segment. The path lengths are those of the ray tracer: `length` and
    [`BeamletOptics.optical_path_length`](@ref) of the beam equal `s_stop` and `opl_stop` of the
    last non-final segments of its branch.

# Examples

```julia
beam = Beam([0, 0, 0], [0, 1, 0])
solve_system!(system, beam)
segs = path_segments(beam; flen = 0.1)
# position of the beam at the geometric path length s along the first branch
s = 0.25
seg = first(filter(sg -> sg.s_start <= s < sg.s_stop, segs))
point = seg.start + (s - seg.s_start) * (seg.stop - seg.start) / (seg.s_stop - seg.s_start)
```
"""
function path_segments(beam::AbstractBeam{T}; flen::Real = 1.0) where {T}
    segments = PathSegment{T}[]
    _path_segments!(segments, beam, 1, T(flen))
    return segments
end

function path_segments(group::AbstractBeamGroup; flen::Real = 1.0)
    bms = beams(group)
    isempty(bms) && return PathSegment{Float64}[]
    T = eltype(position(first(bms)))
    segments = PathSegment{T}[]
    for (i, b) in enumerate(bms)
        _path_segments!(segments, b, i, T(flen))
    end
    return segments
end

"""
    _path_segments!(segments, beam, index, flen, parent = 0, depth = 1, s = 0, opl = 0)

Appends the segments of `beam` and, recursively, of its children to `segments`, see
[`path_segments`](@ref). `parent` is the index of the segment the beam continues from, and `s`,
`opl` are the geometric and optical path length accumulated there.
"""
function _path_segments!(segments::Vector{PathSegment{T}}, beam::AbstractBeam, index::Int, flen::T,
        parent::Int = 0, depth::Int = 1, s::T = zero(T), opl::T = zero(T)) where {T}
    for ray in rays(_chief_beam(beam))
        isect = intersection(ray)
        final = isnothing(isect)
        len = final ? flen : T(length(isect))
        p0 = Point3{T}(position(ray))
        n = refractive_index(ray)
        push!(segments,
            (start = p0, stop = p0 + len * Point3{T}(direction(ray)),
                s_start = s, s_stop = s + len, opl_start = opl, opl_stop = opl + n * len,
                λ = T(wavelength(ray)), depth, parent, beam = index, final))
        parent = length(segments)
        s += len
        opl += n * len
    end
    for child in children(beam)
        _path_segments!(segments, child, index, flen, parent, depth + 1, s, opl)
    end
    return segments
end
