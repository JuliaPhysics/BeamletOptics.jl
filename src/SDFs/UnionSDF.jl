"""
    UnionSDF{T, TT <: Tuple} <: AbstractSDF{T}

This SDF represents the merging of two or more SDFs. If the constituent SDFs do not overlap
(they can and should touch) the resulting SDF should be still exact if the constituent SDFs
are exact.

The intended way to construct these is not explicitely but by just adding two `AbstractSDFs`
using the regular `+` operator.

```julia
s1 = SphereSDF(1.0)
translate3d!(s1, Point3(0, 1.0, 0.0))

s2 = SphereSDF(1.0)

# will result in a SDF with two spheres touching each other.
s_merged = s1 + s2
```

"""
mutable struct UnionSDF{T <: Number, TT <: Tuple{Vararg{AbstractSDF{T}}}} <: AbstractCompositeSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    sdfs::TT
end

operands(u::UnionSDF) = u.sdfs

"""
    thickness(union)

Calculates the thickness of a union of [`AbstractLensSDF`](@ref)s.
"""
function thickness(u::UnionSDF{T}) where T
    t = zero(T)
    for s in u.sdfs
        if hasmethod(thickness, Tuple{typeof(s)})
            t += thickness(s)
        end
    end
    return t
end

function UnionSDF{T}(sdfs::Vararg{AbstractSDF{T}, N}) where {T, N}
    _check_kinematic_members(sdfs)
    UnionSDF{T, typeof(sdfs)}(
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I),
        Point3{T}(zero(T)),
        sdfs
        )
end

"""
    bounding_sphere_of(u::UnionSDF)

Returns the sphere that encloses the bounding spheres of all operands in their current poses, or
`NoBoundingSphere()` if one of the operands has none.
"""
bounding_sphere_of(u::UnionSDF) = SingleBoundingSphere(MultiBoundingSphere(map(bounding_sphere_of, u.sdfs)))

function sdf(s::UnionSDF, pos)
    # sdf to world transform handled by sub-SDFs
    return minimum(sdf(_sdf, pos) for _sdf in s.sdfs)
end

Base.:+(s1::AbstractSDF{T}, s2::AbstractSDF{T}) where T = UnionSDF{T}(s1, s2)
Base.:+(union::UnionSDF{T}, sdf::AbstractSDF{T}) where T = UnionSDF{T}(union.sdfs..., sdf)
Base.:+(sdf::AbstractSDF{T}, union::UnionSDF{T}) where T = UnionSDF{T}(sdf, union.sdfs...)
Base.:+(u1::UnionSDF{T}, u2::UnionSDF{T}) where T = UnionSDF{T}(u1.sdfs..., u2.sdfs...)

# Without this function it is not possible for SDFs encapsulated in a UnionSDF
# to specialize normal3d as always the generic normal3d function is called.
function normal3d(s::UnionSDF, pos)
    ds = map(_sdf -> sdf(_sdf, pos), s.sdfs)
    return normal3d(s.sdfs[_closest_operand(s.sdfs, ds, pos)], pos)
end

"""
    _closest_operand(operands, ds, pos)

Returns the index of the operand closest to `pos`, for the SDF values `ds` of the `operands` at
`pos`. Operands that are equally close up to rounding touch each other tangentially near `pos`,
e.g. the vertex of a concave lens surface touches the face of the lens body. Both give the same
normal, but only the outer one, i.e. the boundary of the union next to `pos`, gives its change,
e.g. its derivative w.r.t. the pose of the lens in automatic differentiation. The outer one is the
closest operand at four points beside `pos`, `ρ` along the tangent plane and `ρ` outward, i.e.
outside of the union also if `pos` lies on its surface.

!!! note
    The offset `ρ` is absolute, 1 µm by default. The tie-break is valid for features that are
    large against `ρ`, i.e. for optics on the mm scale. For features on the µm scale, e.g.
    micro-optics, the probes can land beyond the feature and pick the inner operand. The traced
    values are the same for either pick, as the tied operands give the same normal, but the
    derivatives in automatic differentiation are those of the inner operand then.
"""
function _closest_operand(operands, ds, pos; ρ = 1e-6)
    d_min = minimum(ds)
    tol = Config.get_sdf_raymarch_eps() / 100
    # the usual case, one closest operand, without allocations
    count(d -> d ≤ d_min + tol, ds) == 1 && return findfirst(d -> d ≤ d_min + tol, ds)
    tied = findall(d -> d ≤ d_min + tol, collect(ds))
    n = normal3d(operands[first(tied)], pos)
    t1 = normal3d(Point3(n))
    t2 = cross(n, t1)
    beside(o) = sum(sdf(o, pos + ρ * (t + n)) for t in (t1, -t1, t2, -t2))
    return tied[argmin(map(i -> beside(operands[i]), tied))]
end