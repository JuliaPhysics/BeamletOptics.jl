"""
    UnionSDF{T, TT <: Tuple} <: AbstractSDF{T}

This SDF represents the merging of two or more SDFs. If the constituent SDFs do not overlap
(they can and should touch) the resulting SDF should be still exact if the constituent SDFs
are exact.

The intended way to construct these is not explicitely but by just adding two `AbstractSDFs`
using the regular `+` operator.

```@example
s1 = SphereSDF(1.0)
translate3d!(s1, Point3(0, 1.0, 0.0))

s2 = SphereSDF(1.0)

# will result in a SDF with two spheres touching each other.
s_merged = s1 + s2
```

"""
mutable struct UnionSDF{T <: Number, TT <: Tuple{Vararg{AbstractSDF{T}}}} <: AbstractSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    sdfs::TT
end

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
    UnionSDF{T, typeof(sdfs)}(
        SMatrix{3,3}(one(T)*I),
        SMatrix{3,3}(one(T)*I),
        Point3{T}(zero(T)),
        sdfs
        )
end

function sdf(s::UnionSDF, pos)
    # sdf to world transform handled by sub-SDFs
    return minimum(sdf(_sdf, pos) for _sdf in s.sdfs)
end

Base.:+(s1::AbstractSDF{T}, s2::AbstractSDF{T}) where T = UnionSDF{T}(s1, s2)
Base.:+(union::UnionSDF{T}, sdf::AbstractSDF{T}) where T = UnionSDF{T}(union.sdfs..., sdf)
Base.:+(sdf::AbstractSDF{T}, union::UnionSDF{T}) where T = UnionSDF{T}(sdf, union.sdfs...)
Base.:+(u1::UnionSDF{T}, u2::UnionSDF{T}) where T = UnionSDF{T}(u1.sdfs..., u2.sdfs...)

function translate3d!(u::UnionSDF, offset)
    position!(u, position(u) .+ offset)
    translate3d!.(u.sdfs, Ref(offset))
    return nothing
end

function rotate3d!(u::UnionSDF, axis, θ)
    R = rotate3d(axis, θ)
    # Update group orientation
    orientation!(u, R * orientation(u))
    # Rotate all sub-SDFs around union center
    for s in u.sdfs
        rotate3d!(s, axis, θ)
        v = position(s) - position(u)
        # Translate group around pivot point
        v = (R * v) - v
        translate3d!(s, v)
    end
    return nothing
end

function render_object!(axis, s::UnionSDF)
    for sdf in s.sdfs
        render_object!(axis, sdf)
    end
end
