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
    normal3d(target, reference)

Returns a vector with unit length that is perpendicular to the target and an additional
reference vector. Vector orientation is determined according to right-hand rule.
"""
function normal3d(target, reference)
    n = cross(target, reference)
    return normalize(n)
end

"""
    normal3d(input)

Returns a **random** vector with unit length that is perpendicular to the `input` vector.
"""
function normal3d(input::AbstractArray)
    # Gram-Schmidt method with random init
    new = rand(3)
    # Account for non-normed input vector
    new -= dot(new, input) * input / norm(input)^2
    return normalize(new)
end

function normal3d(input::Point3)
    new = rand(Point3)
    new -= dot(new, input) * input / norm(input)^2
    return normalize(new)
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
        return @SArray [-one(T) zero(T) zero(T);
                        zero(T) -one(T) zero(T);
                        zero(T) zero(T) one(T)]
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
    angle3d(target::AbstractVector, reference::AbstractVector)

Returns the angle between the `target` and `reference` vector in **rad**.
"""
function angle3d(target::AbstractArray{T}, reference::AbstractArray{R}) where {T,R}
    G = promote_type(T,R)
    arg = clamp(dot(target, reference) / (norm(target) * norm(reference)), -one(G), one(G))
    angle = acos(arg)
    return angle
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
function line_plane_distance3d(plane_position::AbstractArray, plane_normal::AbstractArray, line_position::AbstractArray, line_direction::AbstractArray)
    denom = dot(plane_normal, line_direction)
    if abs(denom) > 1e-6
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
    reflection3d(dir, normal)

Calculates the reflection between an input vector `dir` and surface normal vector `normal` in 3D-space
"""
function reflection3d(dir, normal)
    return dir - 2 * dot(dir, normal) * normal
end

"""
    refraction3d(dir, normal, n1, n2)

Calculates the refraction between an input vector `dir` and surface `normal` vector  in 3D-space.
`n1` is the "outside" refractive index and `n2` is the "inside" refractive index.
The function returns the new direction of propagation and a boolean flag to indicate if internal refraction has occured.

Vectors `dir` and `normal` must have **unit length**!

# Total internal reflection

If the critical angle for n1, n2 and the incident angle is reached, the ray is reflected internally instead!

# Arguments

- `dir`: direction vector of incoming ray
- `normal`: surface normal at point of intersection
- `n1`: index of ref. before refraction
- `n2`: index of ref. after refraction
"""
function refraction3d(dir::AbstractArray, normal::AbstractArray, n1::Real, n2::Real)
    # dir and normal must have unit length!
    isapprox(norm(dir), 1) || throw(ArgumentError("dir must have  unit length"))
    isapprox(norm(normal), 1) || throw(ArgumentError("norm must have  unit length"))
    n = n1 / n2
    cosθi = -dot(normal, dir)
    sinθt² = n^2 * (1 - cosθi^2)
    # Check for total reflection
    if sinθt² > 1.0
        return (reflection3d(dir, normal), true)
    end
    cosθt = sqrt(1 - sinθt²)
    return (@. n * dir + (n * cosθi - cosθt) * normal, false)
end

"""
    lensmakers_eq(R1, R2, n)

Calculates the thin lens focal length based on the radius of curvature `R1`/`R2` and the lens refractive index `n`.
If center of sphere is on left then R < 0. If center of sphere is on right then R > 0.
"""
lensmakers_eq(R1, R2, n) = 1 / ((n - 1) * (1 / R1 - 1 / R2))

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

Computes the analytical complex electric field distribution of a stigmatic TEM₀₀ Gaussian beam which is described by:

```math
E(r,z) = {E_0}\\frac{{{w_0}}}{{w(z)}}\\exp\\left( { - \\frac{{{r^2}}}{{w{{(z)}^2}}}} \\right)\\exp\\left(i\\left[ {kz + \\psi + \\frac{{k{r^2}}}{2 R(z)}} \\right] \\right)
```

# Arguments

- `r`: radial distance from beam origin
- `z`: axial distance from beam origin
- `E0`: peak electric field amplitude
- `w0`: waist radius
- `w`: local beam radius
- `k`: wave number, equal to `2π/λ`
- `ψ`: Gouy phase shift (defined as ``-\\text{atan}\\left(\\frac{z}{z_r}\\right)`` !)
- `R`: wavefront curvature, i.e. 1/r (radius of curvature)
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

"Calculates the intensity in [W/m²] for a given complex electric field phasor `E`. Vacuum wave impedance is assumed."
intensity(E::Number, Z=Z_vacuum) = abs2(E) / (2 * Z)

"""
fresnel_coefficients(θ, n)

Calculates the complex Fresnel coefficients for reflection and transmission based on the incident angle `θ` in [rad]
and the refractive index ratio `n = n₂ / n₁`. Returns rₛ, rₚ, tₛ and tₚ.

# Signs

!!! info
    The signs of rₛ, rₚ are based on the definition by Fowles (1975, 2nd Ed. p. 44) and Peatross (2015, 2023 Ed. p. 78)
"""
function fresnel_coefficients(θ::T, n::Number) where T
    C = Complex{T}
    cost = cos(θ)
    n2s2 = sqrt(C(n^2 - sin(θ)^2))
    # Calculate coefficients for reflection/transmission
    rs = (cost - n2s2) / (cost + n2s2)
    rp = (-n^2 * cost + n2s2) / (n^2 * cost + n2s2)
    ts = rs + 1
    tp = 2*n*cost / (n^2 * cost + n2s2)
    return rs, rp, ts, tp
end

function fresnel_coefficients(θ::AbstractArray{T}, n::Number) where T
    rs = Vector{Complex{T}}(undef, length(θ))
    rp = Vector{Complex{T}}(undef, length(θ))
    ts = Vector{Complex{T}}(undef, length(θ))
    tp = Vector{Complex{T}}(undef, length(θ))
    for (i, theta) in enumerate(θ)
        rs[i], rp[i], ts[i], tp[i] = fresnel_coefficients(theta, n)
    end
    return rs, rp, ts, tp
end

is_internally_reflected(rp::Number, rs::Number) = isapprox(abs2(rs), 1, atol=1e-6) && isapprox(abs2(rp), 1, atol=1e-6)
