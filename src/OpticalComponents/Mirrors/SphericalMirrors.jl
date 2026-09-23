"""
    SphericalMirror(radius, thickness, diameter; hole_diameter=nothing)

Constructs a concave spherical [`Mirror`](@ref) with perfect reflectivity (R = 1).
The reflecting surface is modeled as a [`BeamletOptics.UnionSDF`](@ref) of a concave spherical surface
and a plano substrate, combining [`BeamletOptics.ConcaveSphericalSurfaceSDF`](@ref) and [`BeamletOptics.PlanoSurfaceSDF`](@ref).

# Inputs

- `radius`: spherical surface radius of curvature \\[m\\]
- `thickness`: substrate thickness \\[m\\]
- `diameter`: mirror outer diameter \\[m\\]
- `hole_diameter`: diameter of the central through-hole \\[m\\], no hole if `nothing` (default)
"""
function SphericalMirror(radius::Real, thickness::Real, diameter::Real; hole_diameter::Union{Real, Nothing} = nothing)
    cylinder = PlanoSurfaceSDF(thickness, diameter)
    concave = ConcaveSphericalSurfaceSDF(abs(radius), diameter)
    shape = concave + cylinder
    if hole_diameter === nothing
        return Mirror(shape)
    end
    T = float(promote_type(typeof(radius), typeof(thickness), typeof(diameter)))
    r_sph = T(abs(radius))
    d_half = T(diameter) / 2
    sag_max = r_sph >= d_half ? r_sph - sqrt(r_sph^2 - d_half^2) : d_half
    return Mirror(_pierce_substrate(shape, diameter, -sag_max, T(thickness), hole_diameter))
end

# Former name, kept for backwards compatibility
Base.@deprecate_binding ConcaveSphericalMirror SphericalMirror
