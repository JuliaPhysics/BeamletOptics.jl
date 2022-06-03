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
    rotate3d(reference::Ray)

Returns the rotation matrix that will rotate a vector around the reference axis at an angle
θ in radians. Vector length is maintained. Rotation in clockwise direction?
"""
function rotate3d(reference::Vector, θ)
    cost = cos(θ)
    sint = sin(θ)
    ux, uy, uz = reference
    R = [
        cost+ux^2*(1-cost) ux*uy*(1-cost)-uz*sint ux*uz*(1-cost)+uy*sint;
        uy*ux*(1-cost)+uz*sint cost+uy^2*(1-cost) uy*uz*(1-cost)-ux*sint;
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
        return Matrix(1.0I,3,3)
    end
    cosA = dot(start, target)
    k = 1/(1+cosA)
    R = [
        rx^2*k+cosA rx*ry*k+rz rx*rz*k-ry;
        ry*rx*k-rz ry^2*k+cosA ry*rz*k+rx;
        rz*rx*k+ry rz*ry*k-rx rz^2*k+cosA
    ]
    return R
end