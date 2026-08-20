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

function bounding_sphere(u::UnionSDF{T}) where T
    c_merged = Point3{T}(0)
    r_merged = zero(T)
    initialized = false
    for s in u.sdfs
        bs = bounding_sphere(s)
        if bs === nothing
            return nothing
        end
        c_local, r = bs
        c_world = orientation(s) * c_local + position(s)
        if !initialized
            c_merged = c_world
            r_merged = r
            initialized = true
        else
            c2, r2 = c_world, r
            d = norm(c_merged - c2)
            if d + r2 <= r_merged
                continue
            elseif d + r_merged <= r2
                c_merged = c2
                r_merged = r2
            else
                new_r = (d + r_merged + r2) / 2
                if d > 0
                    c_merged = c_merged + (new_r - r_merged) * ((c2 - c_merged) / d)
                end
                r_merged = new_r
            end
        end
    end
    if !initialized
        return nothing
    end
    c_local_union = _world_to_sdf(u, c_merged)
    return (c_local_union, r_merged)
end

# Conversion helper for fields
convert_field(::Type{T}, x::Number) where T = T(x)
convert_field(::Type{T}, x::Point3) where T = Point3{T}(x)
convert_field(::Type{T}, x::Point2) where T = Point2{T}(x)
convert_field(::Type{T}, x::SMatrix{N, M, S, L}) where {T, N, M, S, L} = SMatrix{N, M, T, L}(x)
convert_field(::Type{T}, x::Tuple{Vararg{AbstractSDF}}) where T = map(s -> convert(AbstractSDF{T}, s), x)
convert_field(::Type{T}, x) where T = x

# Reconstruct helper for single-parameter types (default fallback)
reconstruct_sdf(::Type{T}, ::Type{SDFType}, args...) where {T, SDFType} = (Base.typename(SDFType).wrapper){T}(args...)

# Specialization for UnionSDF (multi-parameter)
reconstruct_sdf(::Type{T}, ::Type{<:UnionSDF}, dir, transposed_dir, pos, sdfs) where T =
    UnionSDF{T, typeof(sdfs)}(dir, transposed_dir, pos, sdfs)

# Base conversion method
function Base.convert(::Type{AbstractSDF{T}}, s::SDFType) where {T, S, SDFType <: AbstractSDF{S}}
    if T == S
        return s
    end
    new_fields = map(fieldnames(SDFType)) do name
        convert_field(T, getfield(s, name))
    end
    return reconstruct_sdf(T, SDFType, new_fields...)
end

function Base.:+(s1::AbstractSDF{T1}, s2::AbstractSDF{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    return UnionSDF{T}(convert(AbstractSDF{T}, s1), convert(AbstractSDF{T}, s2))
end
function Base.:+(union::UnionSDF{T1}, sdf::AbstractSDF{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    converted_sdfs = map(s -> convert(AbstractSDF{T}, s), union.sdfs)
    return UnionSDF{T}(converted_sdfs..., convert(AbstractSDF{T}, sdf))
end
function Base.:+(sdf::AbstractSDF{T1}, union::UnionSDF{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    converted_sdfs = map(s -> convert(AbstractSDF{T}, s), union.sdfs)
    return UnionSDF{T}(convert(AbstractSDF{T}, sdf), converted_sdfs...)
end
function Base.:+(u1::UnionSDF{T1}, u2::UnionSDF{T2}) where {T1, T2}
    T = promote_type(T1, T2)
    c1 = map(s -> convert(AbstractSDF{T}, s), u1.sdfs)
    c2 = map(s -> convert(AbstractSDF{T}, s), u2.sdfs)
    return UnionSDF{T}(c1..., c2...)
end

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

# Without this function it is not possible for SDFs encapsulated in a UnionSDF
# to specialize normal3d as always the generic normal3d function is called.
function normal3d(s::UnionSDF, pos)
    # find the closes sub-sdf and call its normal method
    idx = argmin(sdf(_sdf, pos) for _sdf in s.sdfs)

    return normal3d(s.sdfs[idx], pos)
end

function surface_tag(u::UnionSDF, point)
    p_local = _world_to_sdf(u, point)
    n_local = transposed_orientation(u) * normal3d(u, point)
    return surface_tag(u, p_local, n_local)
end

function surface_tag(u::UnionSDF, local_p, local_n)
    # FIXME face_id not defined
    return face_id(u, local_n)
end