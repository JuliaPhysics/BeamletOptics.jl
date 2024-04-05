"""
    PolarizedRay{T} <: AbstractRay{T}

A ray type to model the propagation of an electric field vector based on the publication:

**Yun, Garam, Karlton Crabtree, and Russell A. Chipman. "Three-dimensional polarization ray-tracing calculus I: definition and diattenuation." Applied Optics 50.18 (2011): 2855-2865.**

The geometrical ray description is identical to the standard [`Ray`](@ref). The polarization interaction can be described in local s-p-coordinates 
but must be transformed into global coordinates using the method described in the publication above, see also [`_calculate_global_E0`](@ref).

# Fields

- `id`: a UUID4 that uniquely identifies the `Ray`
- `pos`: a point in R³ that describes the `Ray` origin
- `dir`: a normalized vector in R³ that describes the `Ray` direction
- `intersection`: refer to [`Intersection`](@ref)
- `λ`: wavelength in [m]
- `n`: refractive index along the beam path
- `E0`: complex-valued 3-tuple to represent the electric field in global coordinates

# Jones matrices

In local coordinates the Jones matrices in the case of reflection/refraction are defined as 

- reflection: [-rₛ 0; 0 rₚ]
- transmission: [tₛ 0; 0 tₚ]

where r and t are the complex-valued Fresnel coefficients (see also [`fresnel_coefficients`](@ref)).

# Additional information

!!! warning "Field vector"
    It is assumed that the electric field vector ``E_0`` stays orthogonal to the direction of propagation throughout the optical system.

!!! warning "Intensity"
    E0 can not be converted into an [`intensity`](@ref) value, since a single `PolarizedRay` can not directly model the change in intensity during imaging by an optical system.
"""
mutable struct PolarizedRay{T} <: AbstractRay{T}
    id::UUID
    pos::Point3{T}
    dir::Point3{T}
    intersection::Nullable{Intersection{T}}
    λ::T
    n::T
    E0::Point3{Complex{T}}
    function PolarizedRay{T}(id, pos, dir, intersection, λ, n, E0) where T
        # Crucial check: E0 orthogonal to ray dir.
        if !isorthogonal(dir, E0, atol=1e-14)
            error("Ray direction and field vector E0 must be orthogonal!")
        end
        return new{T}(id, pos, dir, intersection, λ, n, E0)
    end
end

electric_field(ray::PolarizedRay) = ray.E0
electric_field!(ray::PolarizedRay, new) = (ray.E0 = new)

"""
    PolarizedRay(pos, dir, λ = 1000e-9, E0 = [1, 0, 0])

1 V/m in x-dir.
"""
function PolarizedRay(pos::AbstractArray{P},
        dir::AbstractArray{D},
        λ = 1000e-9,
        E0 = [electric_field(1), 0, 0]) where {P <: Real, D <: Real}
    F = promote_type(P, D)
    dir = normalize(dir)
    return PolarizedRay{F}(uuid4(),
        Point3{F}(pos),
        Point3{F}(dir),
        nothing,
        λ,
        F(1),
        E0)
end

"""
    _calculate_global_E0(in_dir, out_dir, J, E0)

Calculates the resulting polarization vector as per the publication by Yun et al. for each surface interaction.
If the `in`- and `out`-directions of propagation are parallel, an arbitrary basis is chosen for the s- and p-components.

# Arguments
- `in_dir`: propagation direction before surface interaction
- `out_dir`: propagation direction after surface interaction
- `J`: Jones matrix extended to 3x3, e.g. [-rₛ 0 0; 0 rₚ 0; 0 0 1] for reflection
- `E0`: Polarization vector before surface interaction
"""
function _calculate_global_E0(in_dir::AbstractArray, out_dir::AbstractArray, J::AbstractArray, E0::AbstractArray)
    s = cross(in_dir, out_dir)
    # Test if in and out dir. are parallel
    if norm(s) ≈ 0
        # Choose arbitrary s-, p-basis.
        s = normal3d(in_dir)
    end
    s = normalize(s)
    # Calculate transforms
    p1 = cross(in_dir, s)
    p2 = cross(out_dir, s)
    O_in = [s'; p1'; in_dir']
    O_out = [s p2 out_dir]
    # Calculate new E0
    return O_out * J * O_in * E0
end