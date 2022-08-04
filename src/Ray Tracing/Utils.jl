"""
    orthogonal3d(target::Vector, reference::Vector)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation is determined according to right-hand rule.
"""
function orthogonal3d(target::Vector, reference::Vector)
    n = cross(target, reference)
    n /= norm(n)
    return n
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
    start /= norm(start)
    target /= norm(target)
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

Returns the angle between the `target` and `reference` vector in **rad**. Also logs the angle in debug mode (in degrees).
"""
function angle3d(target::Vector, reference::Vector)
    angle = acos(dot(target, reference) / (norm(target) * norm(reference)))
    @debug "Angle is $(angle*180/π)°"
    return angle
end

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
    c[1] = a[2] * b[3] - a[3] * b[2]
    c[2] = a[3] * b[1] - a[1] * b[3]
    c[3] = a[1] * b[2] - a[2] * b[1]
    return nothing
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

