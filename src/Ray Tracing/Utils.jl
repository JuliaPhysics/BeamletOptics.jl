"""
    norm3d(v)

Computes the euclidic (p=2) norm of the input `v`ector.
"""
function norm3d(v)
    return @inbounds sqrt(v[1]^2 + v[2]^2 + v[3]^2)
end

"""
    norm3d(v)

Mutation function that changes  the `v`ector length to be ``||v|| = 1``. See also `normalize3d`\\
Note that there is a limit to the equality due to machine precision.
"""
function normalize3d!(v)
    v ./= norm3d(v)
    return nothing
end

"""
    normalize3d(v)

Returns a normalized version of `v`. See also `normalize3d!`.
"""
function normalize3d(v)
    n = copy(v)
    normalize3d!(n)
    return n
end

"""
    orthogonal3d(target::Vector, reference::Vector)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation is determined according to right-hand rule.
"""
function orthogonal3d(target::Vector, reference::Vector)
    n = cross(target, reference)
    return normalize3d(n)
end

"""
    rotate3d(reference::Vector, θ)

Returns the rotation matrix that will rotate a vector around the reference axis at an angle
θ in radians. Vector length is maintained. Rotation in clockwise direction?
"""
function rotate3d(reference::Vector, θ)
    cost = cos(θ)
    sint = sin(θ)
    ux, uy, uz = reference
    R = [
        cost+ux^2*(1-cost) ux*uy*(1-cost)-uz*sint ux*uz*(1-cost)+uy*sint
        uy*ux*(1-cost)+uz*sint cost+uy^2*(1-cost) uy*uz*(1-cost)-ux*sint
        uz*ux*(1-cost)-uy*sint uz*uy*(1-cost)+ux*sint cost+uz^2*(1-cost)
    ]
    return R
end

"""
    align3d(start::Vector, target::Vector)

Returns the rotation matrix R that will align the start vector to be parallel to the target vector.
Based on ['Avoiding Trigonometry'](https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724) by Íñigo Quílez. The resulting matrix
was transposed due to column/row major issues. Vector length is maintained. This function is very fast.
"""
function align3d(start::Vector, target::Vector)
    normalize3d!(start)
    normalize3d!(target)
    rx, ry, rz = cross(target, start)
    # if start and target are already (almost) parallel return unity
    if (abs(rx) < 1e-9) & (abs(ry) < 1e-9) & (abs(rz) < 1e-9)
        return Matrix(1.0I, 3, 3)
    end
    cosA = dot(start, target)
    k = 1 / (1 + cosA)
    R = [
        rx^2*k+cosA rx*ry*k+rz rx*rz*k-ry
        ry*rx*k-rz ry^2*k+cosA ry*rz*k+rx
        rz*rx*k+ry rz*ry*k-rx rz^2*k+cosA
    ]
    return R
end

"""
    angle3d(target::Vector, reference::Vector)

Returns the angle between the `target` and `reference` vector in **rad**.
"""
function angle3d(target::Vector{T}, reference::Vector{T}) where T    
    arg = clamp(fast_dot3d(target, reference) / (norm3d(target) * norm3d(reference)), -one(T), one(T))
    angle = acos(arg)
    return angle
end

angle3d(target::Vector{T}, reference::Vector{V}) where {T, V} = angle3d(promote(target, reference)...)

"""
    fast_dot3d(a, b)

SIMD-accelerated version of the dot product between the vectors `a` and `b`.\\
Assumes a vector dimension of E=3.
"""
function fast_dot3d(a, b)
    c = zero(eltype(a))
    @simd for i = 1:3
        @inbounds c += a[i] * b[i]
    end
    return c
end

"""
    fast_cross3d!(c, a, b)

Specific implementation of the cross product for vector dimensions E=3.\\
Mutates the result in `c` for the cross product of `a` with `b`.
"""
function fast_cross3d!(c, a, b)
    @inbounds c[1] = a[2] * b[3] - a[3] * b[2]
    @inbounds c[2] = a[3] * b[1] - a[1] * b[3]
    @inbounds c[3] = a[1] * b[2] - a[2] * b[1]
    return nothing
end

"""
    fast_cross3d(a, b)

Non-mutating version of fast_cross3d!().
"""
function fast_cross3d(a, b)
    c = zeros(eltype(a), 3)
    fast_cross3d!(c, a, b)
    return c
end

"""
    fast_dot3d(a, b)

SIMD-accelerated version of vector subtraction `b` from `a`.\\
Assumes a vector dimension of E=3 and mutates the data in `c`.
"""
function fast_sub3d!(c, a, b)
    @simd for i = 1:3
        @inbounds c[i] = a[i] - b[i]
    end
    return nothing
end

"""
    line_point_distance3d(pos, dir, point)

Computes the shortes distance between a line described by `pos`+t*`dir` and a `point` in 3D.
This function is slow and should be used only for debugging purposes.
"""
function line_point_distance3d(pos, dir, point)
    return norm3d(SCDI.fast_cross3d(pos .- point, dir)) / norm3d(dir)
end

"""
    reflection3d(dir, normal)

Calculates the reflection between an input vector `dir` and surface normal vector `normal` in 3D-space
"""
function reflection3d(dir, normal)
    return dir .- 2 .* fast_dot3d(dir, normal) .* normal
end

"""
    refraction3d(dir, normal, n1, n2)

Calculates the refraction between an input vector `dir` and surface `normal` vector  in 3D-space.\\
`n1` is the "outside" refractive index and `n2` is the "inside" refractive index.\\
Vectors `dir` and `normal` must have **unit length**!
"""
function refraction3d(dir, normal, n1, n2)
    # dir and normal must have unit length!
    @assert isapprox(norm3d(dir), 1)
    @assert isapprox(norm3d(normal), 1)
    n = n1 / n2
    cosθi = -dot(normal, dir)
    sinθt² = n^2 * (1 - cosθi^2)
    # Check for total reflection
    if sinθt² > 1.0
        return reflection3d(dir, normal)
    end
    cosθt = sqrt(1 - sinθt²)
    return @. n * dir + (n * cosθi - cosθt) * normal
end