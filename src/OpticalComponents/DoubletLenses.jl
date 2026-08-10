abstract type AbstractDoubletRefractiveOptic{
    T,
    F <: AbstractShape{T},
    B <: AbstractShape{T},
    N1 <: RefractiveIndex,
    N2 <: RefractiveIndex
} <: AbstractObjectGroup{T} end

"""
    DoubletLens

Represents a two-component cemented doublet lens with two respective refractive indices `n = n(λ)`.
See also [`SphericalDoubletLens`](@ref).

# Fields

- `front`: front [`Lens`](@ref) component (or [`CoatedRefractive`](@ref))
- `back`: back [`Lens`](@ref) component (or [`CoatedRefractive`](@ref))

# Additional information

!!! warning "Air gap"
    This component type strongly assumes that both lenses are mounted fully flush with respect to each other. 
    Gaps between the components might lead to incorrect results.
"""
struct DoubletLens{T, F <: AbstractRefractiveOptic{T}, B <: AbstractRefractiveOptic{T}} <: AbstractObjectGroup{T}
    front::F
    back::B

    function DoubletLens(front::F, back::B; front_coating = nothing, back_coating = nothing) where {T, F <: AbstractRefractiveOptic{T}, B <: AbstractRefractiveOptic{T}}
        f_coated = !isnothing(front_coating) ? CoatedRefractive(front; front = front_coating) : front
        b_coated = !isnothing(back_coating) ? CoatedRefractive(back; back = back_coating) : back

        # Check if they are flush and aligned along the optical axis
        rel_pos = position(b_coated) - position(f_coated)
        opt_axis = orientation(f_coated)[:, 2]
        expected_offset = thickness(f_coated)
        if !isapprox(rel_pos, expected_offset * opt_axis, atol=1e-5)
            @warn "Doublet components are not aligned flush along the optical axis (expected offset $(expected_offset)m along orientation axis, got $(rel_pos)m). Ensure they were created in the default orientation before assembly. This may cause tracing errors at coincident boundaries."
        end
        if !isapprox(orientation(f_coated), orientation(b_coated), atol=1e-5)
            @warn "Doublet components do not have matching orientation. Ensure they were created in the default orientation before assembly. This may cause tracing errors at coincident boundaries."
        end
        return new{T, typeof(f_coated), typeof(b_coated)}(f_coated, b_coated)
    end
end

shape_trait_of(::DoubletLens) = MultiShape()

shape(dl::DoubletLens) = (dl.front, dl.back)

objects(dl::DoubletLens) = (dl.front, dl.back)

Base.position(dl::DoubletLens) = position(dl.front)
orientation(dl::DoubletLens) = orientation(dl.front)

thickness(dl::DoubletLens) = thickness(shape(dl.front)) + thickness(shape(dl.back))

"""
    SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2; front_coating=nothing, back_coating=nothing)

Generates a two-component "cemented" doublet lens consisting of two spherical lenses.
For radii sign definition, refer to the [`SphericalLens`](@ref) constructor.

# Arguments

- `r1`: radius of curvature for first surface
- `r2`: radius of curvature for second (cemented) surface
- `r3`: radius of curvature for third surface
- `l1`: first lens thickness
- `l2`: second lens thickness
- `d`: lens diameter
- `n1`: first lens [`RefractiveIndex`](@ref)
- `n2`: second lens [`RefractiveIndex`](@ref)
- `front_coating`: optional coating model for front surface
- `back_coating`: optional coating model for back surface
"""
function SphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2; front_coating = nothing, back_coating = nothing)
    # Generate "cemented" front and back spherical lenses
    front = SphericalLens(r1, r2, l1, d, n1)
    back = SphericalLens(r2, r3, l2, d, n2)
    # Move doublet parts into position
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back; front_coating = front_coating, back_coating = back_coating)
end

"""
    AsphericalDoubletLens(r1, r2, r3, l1, l2, d, n1, n2; k1=0, k2=0, k3=0, front_coating=nothing, back_coating=nothing)

Generates a two-component "cemented" doublet lens consisting of two aspherical lens elements with conic constants `k1`, `k2`, and `k3`.
"""
function AsphericalDoubletLens(
        r1, r2, r3, l1, l2, d, n1, n2;
        k1 = 0, k2 = 0, k3 = 0,
        front_coating = nothing, back_coating = nothing
)
    s1 = EvenAsphericalSurface(r1, d, k1, Float64[0.0])
    s2 = EvenAsphericalSurface(r2, d, k2, Float64[0.0])
    s3 = EvenAsphericalSurface(r3, d, k3, Float64[0.0])
    n1_func = (n1 isa Real) ? (λ -> n1) : test_refractive_index_function(n1)
    n2_func = (n2 isa Real) ? (λ -> n2) : test_refractive_index_function(n2)
    front = Lens(s1, s2, l1, n1_func)
    back = Lens(s2, s3, l2, n2_func)
    translate3d!(back, [0, thickness(shape(front)), 0])
    return DoubletLens(front, back; front_coating = front_coating, back_coating = back_coating)
end