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
function normal3d(target::AbstractVector, reference::AbstractVector)
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

Calculates the reflection between an input vector `dir` and surface `normal` vector in R³.
Vectors `dir` and `normal` must have **unit length**!
"""
function reflection3d(dir, normal)
    return dir - 2 * dot(dir, normal) * normal
end

"""
    refraction3d(dir, normal, n1, n2)

Calculates the refraction between an input vector `dir` and surface `normal` vector in R³.
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
    isapprox(norm(normal), 1) || throw(ArgumentError(lazy"norm must have  unit length: $(norm(normal))"))
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


"""
    sag(r::Real, l::Real)

Calculates the sag of a cut circle with radius `r` and chord length `l`
"""
sag(r::Real, l::Real) = r - sqrt(r^2 - 0.25 * l^2)

function check_sag(r, d)
    if abs(2r) < d
        throw(ArgumentError("Radius of curvature (r = $(r)) must be ≥ than half the diameter (d = $(d)) or an illegal shape results!"))
    else
        return nothing
    end
end

## Refractive index utils
"""
    DiscreteRefractiveIndex{T}

Represents a incomplete set of dispersion data where for each exact wavelength one refractive index value is stored in the `data` field.
Can be called like a function `n = n(λ)`. Does not interpolate between data points.
Refer to [`RefractiveIndex`](@ref) for more information.
"""
struct DiscreteRefractiveIndex{T}
    data::Dict{T, T}
end

"""
    DiscreteRefractiveIndex(λs, n)

Creates a [`DiscreteRefractiveIndex`](@ref) dictionary where each wavelength in `λs` is mapped onto an exact exact refractive index in `ns`.

# Inputs

- `λs`: array of wavelengths
- `ns`: array of refractive indices
"""
function DiscreteRefractiveIndex(λs::AbstractArray{L}, ns::AbstractArray{N}) where {L, N}
    if length(λs) != length(ns)
        throw(ArgumentError("Number of wavelengths must match number of ref. indices"))
    end
    T = promote_type(L, N)
    d = Dict(zip(λs, ns))
    return DiscreteRefractiveIndex{T}(d)
end

(dri::DiscreteRefractiveIndex)(λ) = dri.data[λ]

"[`DiscreteRefractiveIndex`](@ref) passes test by default"
test_refractive_index_function(::DiscreteRefractiveIndex) = nothing

"""
    test_refractive_index_function(input)

Tests if `input` is callable with a single `Real` argument for the wavelength `λ` and
returns a single `Real` value for the refractive index `n`.
"""
function test_refractive_index_function(input)
    # Test function compat for the following types of λ
    Ts = [Int, Float32, Float64]
    try
        for T in Ts
            # Test if input accepts types
            answer = input(one(T))
            # Test if input returns single real value
            if !isa(answer, Real)
                error()
            end
        end
    catch
        error_msg = "Ref. index must be callable with a single Real argument and return a single real result."
        throw(ArgumentError(error_msg))
    end
    return nothing
end

"""
    RefractiveIndex

Union type that represents valid means to pass a refractive index `n` to e.g. [`AbstractObject`](@ref)s.
The core assumption is that:

1. the refractive index is callable with a **single** `Number` argument `λ` to represent the wavelength in [m]
2. the return value is a **single** `Number` value for the refractive index

Refer to e.g. [`DiscreteRefractiveIndex`](@ref).
"""
const RefractiveIndex = Union{Function, DiscreteRefractiveIndex}

"""
    list_subtypes(T::Type; max_depth::Int=5)

Prints a tree of all subtypes, e.g. `list_subtypes(AbstractObject)`.
Maximum exploration depth can be limited by passing `max_depth`.
Returns the total number of types encountered.
"""
function list_subtypes(parent::Type, is_last=true, prefix="", depth=1; max_depth=10)
    # Determine children of parent type
    children = subtypes(parent)
    num_children = length(children)
    # Print out parent type information
    connector = is_last ? "└── " : "├── "
    msg = prefix * connector * string(parent)
    println(msg)
    # Return to parent function if childless
    if iszero(num_children)
        return 0
    end
    # Extend prefix string based on layer depth
    if is_last
        prefix = prefix * "    "
    else
        prefix = prefix * "│   "
    end
    # Check if max tree depth has been reached
    if depth > max_depth
        msg = prefix * "└── " * "$num_children more ..."
        println(msg)
        return num_children
    end
    # Delve deeper into tree
    for child in children[1:end-1]
        num_children += list_subtypes(child, false, prefix, depth+1; max_depth)
    end
    num_children += list_subtypes(last(children), true, prefix, depth+1; max_depth)
    # If tree/branch has been resolved, return total number of types found.
    if depth == 1
        msg = "\n" * "At least $num_children types have been found."
        println(msg)
    end
    return num_children
end

"""
    countlines_in_dir(dir)

Counts the number of lines of all `.jl` files in `dir`.
"""
function countlines_in_dir(dir::String)
    l = 0
    items = readdir(dir)
    for item in items
        path = joinpath(dir, item)
        if isfile(path) && endswith(path, ".jl")
            li = countlines(path)
            l += li
            println("$li lines of code in $path")
        end
        if isdir(path)
            l += countlines_in_dir(path)
        end
    end
    return l
end

"""
    ToDO
"""
function find_zero_bisection(f, a, b; tol=1e-10, max_iter=1000)
    fa = f(a)
    fb = f(b)
    if sign(fa) == sign(fb)
        error("Bisection requires a sign change: f(a)=$(fa), f(b)=$(fb)")
    end
    for i in 1:max_iter
        mid = (a + b) / 2
        fmid = f(mid)
        if abs(fmid) < tol
            return mid
        end
        if sign(fa) == sign(fmid)
            a = mid
            fa = fmid
        else
            b = mid
            fb = fmid
        end
    end
    error("Bisection did not converge after $max_iter iterations")
end

"""
    numerical_aperture(θ, n=1)

Returns the `NA` for a opening half-angle `θ` and scalar ref. index `n`.
For more information refer to [this website](https://www.rp-photonics.com/numerical_aperture.html).
"""
numerical_aperture(θ::Real, n::Real=1.0) = n*sin(θ)