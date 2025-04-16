"""
    Photodetector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}

Represents a **flat** rectangular or quadratic surface in R³ that is the active surface of a photodetector.
The active surface is discretized in the local R² x-y-coordinate system.
Field contributions Eᵢ are added by the corresponding [`interact3d`](@ref) method.

# Fields

- `shape`: geometry of the active surface, must represent 2D-`field` in `x` any `y` dimensions
- `x`: linear range of local x-coordinates
- `y`: linear range of local y-coordinates
- `field`: `size(x)` by `size(y)` matrix of complex values to store superposition E₀

# Additional information

!!! warning "Reset behavior"
    The `Photodetector` must be reset between each call of [`solve_system!`](@ref) in order to
    overwrite previous results using the [`reset_detector!`](@ref) function.
    Otherwise, the current result will be added onto the previous result.

!!! info "Supported beams"
    Currently, only the [`GaussianBeamlet`](@ref) is supported.
"""
mutable struct Photodetector{T, S <: AbstractShape{T}} <: AbstractDetector{T, S}
    const shape::S
    x::LinRange{T, Int}
    y::LinRange{T, Int}
    field::Matrix{Complex{T}}
end

"""
    Photodetector(width, n)

Spawns a quadratic rectangular 2D [`Photodetector`](@ref) that is aligned with the **positive y-axis**.
Refer to the type docs for more information.

# Inputs:

- `width`: edge length in [m]
- `n`: field discretization factor, higher results in more computational cost
"""
function Photodetector(width::T, n::Int) where {T<:Real}
    shape = QuadraticFlatMesh(width)
    sz = maximum(vertices(shape))
    x = y = LinRange(-sz, sz, n)
    field = zeros(Complex{T}, n, n)
    return Photodetector{T, typeof(shape)}(shape, x, y, field)
end

function interact3d(::AbstractSystem, ::Photodetector, ::B, ::Ray) where {B <: AbstractBeam}
    @warn "Photodetection for $B not implemented"
    return nothing
end

"""
    interact3d(::AbstractSystem, pd::Photodetector, gauss::GaussianBeamlet, ray_id::Int)

Implements the [`Photodetector`](@ref) interaction with a [`GaussianBeamlet`](@ref).
On hit, the scalar E-field of the `gauss` is added to the current PD field matrix.
Tilt and tip between beam and PD surface are considered via projection factors.
"""
function interact3d(
        ::AbstractSystem, pd::Photodetector, gauss::GaussianBeamlet, ray_id::Int)
    # Select final ray of chief beam
    ray = gauss.chief.rays[ray_id]
    # Subtract ray length from optical path l0 (calculated seperately with projection l1)
    l0 = optical_path_length(gauss) - optical_path_length(ray)
    p0 = position(ray)
    d0 = direction(ray)
    # Preallocate transforms
    T = transpose(orientation(shape(pd)))
    p = position(shape(pd))
    # E-field projection scalar reduction factor
    ray_int = intersection(ray)
    isnothing(ray_int) && return nothing

    proj = abs(dot(d0, normal3d(ray_int)))

    # Add current E-field contribution
    Threads.@threads for j in eachindex(pd.y) # FIXME row column major order?
        y = pd.y[j]
        @inbounds for i in eachindex(pd.x)
            x = pd.x[i]
            # Transform point p on PD into world coordinates
            p1 = Point3(
                T[1, 1] * x + T[1, 3] * y + p[1],
                T[2, 1] * x + T[2, 3] * y + p[2],
                T[3, 1] * x + T[3, 3] * y + p[3]
            )
            # Find projection of p1 onto Gaussian optical axis, i.e. local r and z
            l1 = dot(p1 - p0, d0)
            p2 = p0 + l1 * d0
            r = norm(p1 - p2)
            z = l0 + l1
            # Add field contribution, projection factor accounts for beam spot stretching
            pd.field[i, j] += electric_field(gauss, r, z) * sqrt(proj)
        end
    end
    return nothing
end

"""
    interact3d(::AbstractSystem, pd::Photodetector, beam::Beam{T, R}, ray_id::R) where {T, R <: AbstractRay}

Implements the [`Photodetector`](@ref) interaction with a [`Beam`](@ref).
On hit, the scalar E-field of a plane wave is added to the current PD field matrix.
Tilt and tip between beam and PD surface are considered via projection factors.
"""
function interact3d(
    ::AbstractSystem, pd::Photodetector, beam::Beam{T, R}, ray::R) where {T, R <: AbstractRay}

    # Optical path length difference
    l0 = optical_path_length(beam) - optical_path_length(ray)
    p0 = position(ray)
    d0 = direction(ray)

    orient = orientation(shape(pd))
    e1 = @view orient[:, 1]   # local x-axis
    e2 = @view orient[:, 3]   # local y-axis
    origin_pd = position(shape(pd))

    k = 2π / wavelength(ray)

    # Loop over the detector grid points
    Threads.@threads for j in eachindex(pd.y)
        y = pd.y[j]
        @inbounds for i in eachindex(pd.x)
            x = pd.x[i]
            # Convert detector local coordinates (x,y) to world coordinates
            p1 = origin_pd + x * e1 + y * e2

            # Project the vector from the ray's position to the grid point onto d0.
            l1 = dot(p1 - p0, d0)
            # Compute the effective optical path length (z) relative to the flat phase reference.
            z = l0 + l1
            phase = exp(im * k * z)

            pd.field[i, j] += phase
        end
    end
    return nothing
end

intensity(pd::Photodetector) = intensity.(pd.field)

"""
    optical_power(pd::Photodetector)

Calculates the total optical power on `pd` in [W] by integration over the local intensity.
"""
optical_power(pd::Photodetector) = trapz((pd.x, pd.y), intensity(pd))

"""Resets the values currently stored in `pd.field` to zero"""
reset_detector!(pd::Photodetector{T, S}) where {T, S} = (pd.field .= zero(Complex{T}))

"""
    photodetector_resolution!(pd::Photodetector, n::Int)

Sets the resolution of `pd` to `n` × `n`. Note that this resets the current `pd.field`.
"""
function photodetector_resolution!(pd::Photodetector{T, S}, n::Int) where {T, S}
    pd.x = LinRange(pd.x.start, pd.x.stop, n)
    pd.y = LinRange(pd.y.start, pd.y.stop, n)
    pd.field = zeros(Complex{T}, n, n)
    return nothing
end
