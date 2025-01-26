"""
    PlanoConvexAsphericalLens(radius, conic_constant, even_coefficients, d, t, n, md = d)

Creates a plano-aspherical [`Lens`](@ref) with convex shape based on:

- `radius`: front radius
- `conic_constant`: Conic constant of the aspherical surface
- `even_coefficients`: A vector of the even aspheric coefficients
- `d`: lens diameter
- `t`: lens thickness
- `n`: [`RefractiveIndex`](@ref) as a function of λ

!!! note
    Aspheric lenses are somewhat experimental at the moment. Use this feature with some caution.
    Future versions of this package will offer a convenience constructor for abitrary mixed
    lenses, e.g. bi-aspheres, aspheric-spheric lenses, etc. as well as odd aspheres and extended
    aspheres.
"""
function PlanoConvexAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Real)
    return PlanoConvexAsphericalLens(radius, conic_constant, even_coefficients, d, t, x -> n)
end

function PlanoConvexAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::RefractiveIndex)
    test_refractive_index_function(n)
    shape = PlanoConvexAsphericalLensSDF(radius, t, d, conic_constant, even_coefficients)
    return Lens(shape, n)
end

"""
    PlanoConcaveAsphericalLens(radius, conic_constant, even_coefficients, d, t, n, md = d)

Creates a plano-aspherical [`Lens`](@ref) with convex shape based on:

- `radius`: front radius
- `conic_constant`: Conic constant of the aspherical surface
- `even_coefficients`: A vector of the even aspheric coefficients
- `d`: lens diameter
- `t`: lens thickness
- `n`: [`RefractiveIndex`](@ref) as a function of λ
- `md`: mechanical diameter of the lens (defaults to `d`). If this is set to a value `md` > `d`
        an outer flat section will be added to the lens. This can be used to model more realistic
        lenses where this flat section is present for mounting purposes but is also an active
        optical region for extreme lenses (HUDs, AR-wearables, etc.)

!!! note
    Aspheric lenses are somewhat experimental at the moment. Use this feature with some caution.
    Future versions of this package will offer a convenience constructor for abitrary mixed
    lenses, e.g. bi-aspheres, aspheric-spheric lenses, etc. as well as odd aspheres and extended
    aspheres.
"""
function PlanoConcaveAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::Real, md::Real = d)
    return PlanoConcaveAsphericalLens(radius, conic_constant, even_coefficients, d, t, x -> n, md)
end

function PlanoConcaveAsphericalLens(radius::Real, conic_constant::Real, even_coefficients::Vector{<:Real}, d::Real, t::Real, n::RefractiveIndex, md::Real = d)
    test_refractive_index_function(n)
    shape = PlanoConcaveAsphericalLensSDF(radius, t, d, conic_constant, even_coefficients, md)
    return Lens(shape, n)
end