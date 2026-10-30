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

!!! warning
    Where operands touch tangentially, the normal is taken from the operand that is closest at
    probes `ρ` = 1 µm beside the surface, which is valid for features much larger than `ρ`. The
    traced values are unaffected, but on the µm scale the derivatives in automatic
    differentiation may come from the inner operand.
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

Calculates the thickness of a union of [`AbstractLensSDF`](@ref)s as the sum of the thicknesses
of its operands, e.g. 0 for a mount ring. Throws a `MethodError` for an operand without
`thickness`.
"""
function thickness(u::UnionSDF{T}) where T
    t = zero(T)
    for s in u.sdfs
        t += thickness(s)
    end
    return t
end

function UnionSDF{T}(sdfs::Vararg{AbstractSDF{T}, N}) where {T, N}
    check_kinematic_members(sdfs)
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

sdf(s::UnionSDF, pos) = minimum(map(o -> sdf(o, pos), s.sdfs))

Base.:+(s1::AbstractSDF{T}, s2::AbstractSDF{T}) where T = UnionSDF{T}(s1, s2)
Base.:+(union::UnionSDF{T}, sdf::AbstractSDF{T}) where T = UnionSDF{T}(union.sdfs..., sdf)
Base.:+(sdf::AbstractSDF{T}, union::UnionSDF{T}) where T = UnionSDF{T}(sdf, union.sdfs...)
Base.:+(u1::UnionSDF{T}, u2::UnionSDF{T}) where T = UnionSDF{T}(u1.sdfs..., u2.sdfs...)

# Without this function it is not possible for SDFs encapsulated in a UnionSDF
# to specialize normal3d as always the generic normal3d function is called.
function normal3d(s::UnionSDF, pos)
    ds = map(o -> sdf(o, pos), s.sdfs)
    return normal3d(s.sdfs[closest_operand(s, ds, pos)], pos)
end

"""
    closest_operand(u::UnionSDF, ds, pos; ρ = 1e-6)

Returns the index of the operand of `u` closest to `pos`, for its SDF values `ds` at `pos`; ties
go to the operand closest at probes `ρ` \\[m\\] outside.
"""
function closest_operand(u::UnionSDF, ds, pos; ρ = 1e-6)
    operands = u.sdfs
    d_min = minimum(ds)
    tol = Config.get_sdf_raymarch_eps() / 100
    # Single closest operand without allocations
    count(d -> d ≤ d_min + tol, ds) == 1 && return findfirst(d -> d ≤ d_min + tol, ds)
    tied = findall(d -> d ≤ d_min + tol, collect(ds))
    n = normal3d(operands[first(tied)], pos)
    t1 = normal3d(Point3(n))
    t2 = cross(n, t1)
    beside(o) = sum(sdf(o, pos + ρ * (t + n)) for t in (t1, -t1, t2, -t2))
    return tied[argmin(map(i -> beside(operands[i]), tied))]
end