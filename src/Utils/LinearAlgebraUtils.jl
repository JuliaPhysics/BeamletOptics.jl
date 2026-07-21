"""
    isparallel3d(v1, v2)

Tests if `v1` is parallel to `v2`.
"""
function isparallel3d(v1::AbstractArray, v2::AbstractArray)
    return isapprox(abs(dot(normalize(v1), normalize(v2))), 1, atol=eps())
end

"""
    isorthogonal3d(v1, v2; atol=eps())

Tests if `v1` and `v2` are orthogonal. Additional abs. tolerance can be passed via `atol`
"""
function isorthogonal3d(v1::AbstractArray, v2::AbstractArray; atol=eps())
    return isapprox(dot(v1, v2), 0; atol)
end

"""
    normal3d(target, reference)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation is determined according to right-hand rule.
"""
function normal3d(target::AbstractVector, reference::AbstractVector)
    n = cross(target, reference)
    return normalize(n)
end

"""
    normal3d(input)

Returns a vector with unit length that is perpendicular to the `input` vector.
The orientation is chosen deterministically to guarantee reproducible bases.
"""
@inline function _orthogonal_basis_vector(v::SVector{3,T}) where {T}
    if abs(v[1]) > abs(v[2])
        inv_len = inv(sqrt(v[1]^2 + v[3]^2))
        return SVector(-v[3] * inv_len, zero(T), v[1] * inv_len)
    else
        inv_len = inv(sqrt(v[2]^2 + v[3]^2))
        return SVector(zero(T), v[3] * inv_len, -v[2] * inv_len)
    end
end

function normal3d(input::AbstractArray)
    T = float(eltype(input))
    v = normalize(SVector{3,T}(Tuple(input)))
    n = _orthogonal_basis_vector(v)
    # stabilize output type
    if input isa SVector
        return convert(typeof(input), n)
    elseif input isa AbstractVector
        return collect(n)
    else
        return convert(typeof(input), n)
    end
end

function normal3d(input::Point3{T}) where T
    v = normalize(SVector{3,T}(Tuple(input)))
    n = _orthogonal_basis_vector(v)
    return Point3(n...)
end

"""
    rotate3d(reference::Vector, θ)

Returns the rotation matrix that will rotate a vector around the reference axis at an angle
θ in radians. Vector length is maintained. Counter-clockwise rotation in a right-hand coord. system.
"""
function rotate3d(reference::AbstractVector, θ)
    if isnan(θ) || isinf(θ)
        throw(ArgumentError("θ must be real and not Inf or NaN"))
    end
    cost = cos(θ)
    sint = sin(θ)
    ux, uy, uz = reference
    R = @SArray [
        cost+ux^2*(1-cost) ux*uy*(1-cost)-uz*sint ux*uz*(1-cost)+uy*sint
        uy*ux*(1-cost)+uz*sint cost+uy^2*(1-cost) uy*uz*(1-cost)-ux*sint
        uz*ux*(1-cost)-uy*sint uz*uy*(1-cost)+ux*sint cost+uz^2*(1-cost)
    ]
    return R
end

"""
    align3d(start::AbstractVector, target::AbstractVector)

Returns the rotation matrix R that will align the start vector to be parallel to the target vector.
Based on ['Avoiding Trigonometry'](https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724) by Íñigo Quílez. The resulting matrix
was transposed due to column/row major issues. Vector length is maintained. This function is very fast.
"""
function align3d(start::AbstractVector{A}, target::AbstractVector{B}) where {A, B}
    T = promote_type(A, B)
    start = normalize(start)
    target = normalize(target)
    rx, ry, rz = cross(target, start)
    cosA = dot(start, target)
    # if start and target are already (almost) parallel return unity
    if cosA ≈ 1
        return SMatrix{3,3}(one(T)I)
    end
    if cosA ≈ -1
        u = normal3d(start)
        ux, uy, uz = u
        return @SArray [
            2*ux^2-one(T)  2*ux*uy         2*ux*uz
            2*uy*ux        2*uy^2-one(T)  2*uy*uz
            2*uz*ux        2*uz*uy        2*uz^2-one(T)
        ]
    end
    k = 1 / (1 + cosA)
    R = @SArray [
        rx^2*k+cosA rx*ry*k+rz rx*rz*k-ry
        ry*rx*k-rz ry^2*k+cosA ry*rz*k+rx
        rz*rx*k+ry rz*ry*k-rx rz^2*k+cosA
    ]
    return R
end

"""
    angle3d(target::AbstractVector, start::AbstractVector)

Returns the angle between the `target` and `start` vector in **rad**.
"""
function angle3d(target::AbstractArray{T}, start::AbstractArray{R}) where {T,R}
    G = promote_type(T,R)
    arg = clamp(dot(target, start) / (norm(target) * norm(start)), -one(G), one(G))
    angle = acos(arg)
    return angle
end

"""
    angle3d(target::AbstractArray, start::AbstractArray, reference::AbstractArray)

Returns the angle between the `target` and `start` vector in **rad**. In addition, a `reference` axis must be specified.
This axis is used in order to determine the angle sign of rotation according to the **right hand rule**.
"""
function angle3d(target::AbstractArray{T}, start::AbstractArray{S}, reference::AbstractArray{R}) where {T,S,R}
    G = promote_type(T,S,R)
    θ = angle3d(target, start)
    # get angle sign
    c = cross(start, target)
    if dot(c, reference) > zero(G)
        return θ
    else
        return -θ
    end
end

"""
    line_point_distance3d(pos, dir, point)

Computes the shortes distance between a line described by `pos`+t*`dir` and a `point` in 3D.
This function is slow and should be used only for debugging purposes.
"""
function line_point_distance3d(pos, dir, point)
    d = pos - point
    c = cross(d, dir)
    return norm(c) / norm(dir)
end

"""
    line_plane_distance3d(plane_position, plane_normal, line_position, line_direction)

Returns the distance between a line and an infinitely large plane which are characterized by their `position` and `normal`/`direction`.
"""
function line_plane_distance3d(plane_position::AbstractArray, plane_normal::AbstractArray,
        line_position::AbstractArray, line_direction::AbstractArray)
    denom = dot(plane_normal, line_direction)
    if abs(denom) > Config.get_line_plane_intersection_threshold()
        # explicit dot product for perfomance
        c = dot(plane_position - line_position, plane_normal)
        t = c / denom
        return t
    end
    return nothing
end

"""
    isinfrontof(point::AbstractVector, pos::AbstractVector, dir::AbstractVector)

Tests if a `point` is in front of the plane defined by the `pos`ition and `dir`ection vectors.
"""
function isinfrontof(point::AbstractVector, pos::AbstractVector, dir::AbstractVector)
    los = normalize(point - pos)
    if dot(dir, los) ≤ 0
        return false
    else
        return true
    end
end

"""
    base_transform(base, base2=I(3))

Return the base transformation matrix for transforming from vectors given
relative to `base2` into `base`.
"""
base_transform(base, base2=I(3)) = base \ base2
