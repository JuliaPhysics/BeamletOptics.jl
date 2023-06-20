"""
    Nullable{T}

An alias which results in `Union{T, Nothing}` to provide a shorter notation for struct
fields which can containing nothing.
"""
const Nullable{T} = Union{T,Nothing} where {T}

"""
    NullableVector{T}

An alias which results in `Union{Vector{T}, Nothing}` to provide a shorter notation for struct
fields which can containing nothing.
"""
const NullableVector{T} = Union{Vector{T},Nothing} where {T}

"""
    norm3d(v)

Computes the euclidic (p=2) norm of the 3D input `v`ector.
"""
function norm3d(v)
    return @inbounds sqrt(v[1]^2 + v[2]^2 + v[3]^2)
end

"""
    norm3d(v)

Computes the euclidic (p=2) norm of the 2D input `v`ector.
"""
function norm2d(v)
    return @inbounds sqrt(v[1]^2 + v[2]^2)
end

"""
    normalize3d!(v)

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
    normal3d(target, reference)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation is determined according to right-hand rule.
"""
function normal3d(target, reference)
    n = fast_cross3d(target, reference)
    return normalize3d(n)
end

"""
    rotate3d(reference::Vector, θ)

Returns the rotation matrix that will rotate a vector around the reference axis at an angle
θ in radians. Vector length is maintained. Rotation in clockwise direction?
"""
function rotate3d(reference::AbstractVector, θ)
    cost = cos(θ)
    sint = sin(θ)
    ux, uy, uz = reference
    R = @SMatrix [
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
function align3d(start::Vector{A}, target::Vector{B}) where {A, B}
    T = promote_type(A, B)
    normalize3d!(start)
    normalize3d!(target)
    rx, ry, rz = fast_cross3d(target, start)
    cosA = fast_dot3d(start, target)
    # if start and target are already (almost) parallel return unity
    if cosA ≈ 1
        return Matrix{T}(I, 3, 3)
    end
    if cosA ≈ -1
        return [-one(T) zero(T) zero(T); zero(T) -one(T) zero(T); zero(T) zero(T) one(T)]
    end
    k = 1 / (1 + cosA)
    R = @SMatrix [
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
function angle3d(target::Vector{T}, reference::Vector{T}) where {T}
    arg = clamp(fast_dot3d(target, reference) / (norm3d(target) * norm3d(reference)), -one(T), one(T))
    angle = acos(arg)
    return angle
end

angle3d(target::Vector{T}, reference::Vector{V}) where {T,V} = angle3d(promote(target, reference)...)

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
    # WARNING: can cause buffer overflow
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
function line_point_distance3d(pos, dir, point, buf_a=zeros(length(pos)), buf_b=zeros(length(pos)))
    for i in eachindex(pos)
        buf_a[i] = pos[i] - point[i]
    end
    fast_cross3d!(buf_b, buf_a, dir)
    return norm3d(buf_b) / norm3d(dir)
end

"""
    line_plane_distance3d(plane_position, plane_normal, line_position, line_direction)

Returns the distance between a line and an infinitely large plane which are characterized by their `position` and `normal`/`direction`.
"""
function line_plane_distance3d(plane_position::AbstractArray, plane_normal::AbstractArray, line_position::AbstractArray, line_direction::AbstractArray)
    denom = fast_dot3d(plane_normal, line_direction)
    if denom > 1e-6
        # explicit dot product for perfomance
        c = zero(eltype(plane_position))
        pr = line_position
        @simd for i = 1:3
            @inbounds c += (plane_position[i] - pr[i]) * plane_normal[i]
        end
        t = c / denom
        return t
    end
    return nothing
end

"""
    isinfrontof(point::AbstractVector, pos::AbstractVector, dir::AbstractVector)

Tests if a `point` is in front of the plane defined by the `pos`ition and `dir`ection vectors.
"""
function isinfrontof(point::AbstractVector, pos::AbstractVector, dir::AbstractVector, los=similar(pos))
    @inbounds @simd for i in eachindex(pos)
        los[i] = point[i] - pos[i]
    end
    if fast_dot3d(dir, los) <= 0
        return false
    else
        return true
    end
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

"""
    lensmakers_eq(R1, R2, nl, n0=1.0)

Calculates the thin lens focal length based on the radius of curvature `R1`/`R2` and the lens refractive index `nl`.
If center of sphere is on left then R < 0. If center of sphere is on right then R > 0.
"""
lensmakers_eq(R1, R2, nl, n0=1.0) = 1 / ((nl - n0) / n0 * (1 / R1 - 1 / R2))

"""
    base_transform(base, base2=I(3))

Return the base transformation matrix for transforming from vectors given
relative to `base2` into `base`.
"""
base_transform(base, base2=I(3)) = base \ base2

rayleigh_range(λ, w0, M2, n=1) = π * n * w0^2 / λ / M2

beam_waist(z, w0, zr) = w0 * sqrt(1 + (z / zr)^2)

gouy_phase(z, zr) = -atan(z / zr) # some definitions with +

wavefront_curvature(z, zr) = z / (z^2 + zr^2) # 1/r

divergence_angle(λ, w0, M2) = M2 * λ / (π * w0)

wave_number(λ) = 2π / λ

"""
    electric_field(r, z, E0, w0, w, k, ψ, R) -> ComplexF64

Computes the theoretical complex electric field distribution of an unastigmatic Gaussian laser beam.

# Arguments
- `r`: radial distance from beam origin
- `z`: axial distance from beam origin
- `E0`: peak electric field amplitude
- `w0`: waist radius
- `w`: local beam radius
- `k`: wave number, equal to `2π/λ`
- `ψ`: Gouy phase shift (defined as -atan(z/zr) !)
- `R`: wavefront curvature
"""
electric_field(r::Real, z::Real, E0, w0, w, k, ψ, R) = E0 * w0 / w * exp(-r^2 / w^2) * exp(im * (k * z + ψ + (k * r^2 * R) / 2))

function electric_field(r::Real, z::Real, E0, w0, λ, M2=1)
    zr = rayleigh_range(λ, w0, M2)
    w = beam_waist(z, w0, zr)
    k = wave_number(λ)
    ψ = gouy_phase(z, zr)
    R = wavefront_curvature(z, zr)
    return electric_field(r, z, E0, w0, w, k, ψ, R)
end

"""
    electric_field(I::Real, Z=Z_vacuum, ϕ=0)

Calculates the E-field phasor in [V/m] for a given intensity `I` and phase ϕ. Vacuum wave impedance is assumed.
"""
electric_field(I::Real, Z=Z_vacuum, ϕ=0) = sqrt(2 * I * Z) * exp(im * ϕ)

"Calculates the intensity in [W/m²] for a given electric field phasor `E`. Vacuum wave impedance is assumed."
intensity(E::Complex, Z=Z_vacuum) = abs2(E) / (2 * Z)