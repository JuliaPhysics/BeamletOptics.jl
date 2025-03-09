"""
    MeniscusLensSDF

[`AbstractSDF`](@ref)-based representation of a positive or negative meniscus lens.
When constructed, it is assumed that the lens lies at the origin and the optical axis is aligned with the `y`-axis.
Parameters that lead to a sharp lens edge will cause an error.

# Notes

!!! info "Radius of curvature sign convention"
    The ROC is defined to be positive if the center is to the right of the surface. Otherwise it is negative.

# Fields:

- `convex`: the convex part of the lens composite SDF
- `cylinder`: the cylindrical part of the lens composite SDF
- `concave`: the concave part of the lens composite SDF
- `thickness`: lens thickness on the optical axis
"""
mutable struct MeniscusLensSDF{T, S1 <: AbstractLensSDF{T}, S2 <: AbstractLensSDF{T}} <:
               AbstractLensSDF{T}
    dir::SMatrix{3, 3, T, 9}
    transposed_dir::SMatrix{3, 3, T, 9}
    pos::Point3{T}
    convex::S1
    cylinder::PlanoSurfaceSDF{T}
    concave::S2
    thickness::T
end

diameter(ml::MeniscusLensSDF) = diameter(ml.cylinder)
thickness(ml::MeniscusLensSDF) = ml.thickness

function bounding_sphere(ml::MeniscusLensSDF)
    d = diameter(ml.cylinder)
    l = thickness(ml)
    center = Point3(0, l / 2, 0)
    radius = sqrt((l / 2)^2 + (d / 2)^2)
    return (center, radius)
end

function sdf(ml::MeniscusLensSDF, pos)
    p = _world_to_sdf(ml, pos)
    # Add front and mid, subtract back
    return max(min(sdf(ml.convex, p), sdf(ml.cylinder, p)), -sdf(ml.concave, p))
end

function MeniscusLensSDF(r1::R1, r2::R2, l::L, d::D = 1inch) where {R1, R2, L, D}
    T = promote_type(R1, R2, L, D)
    # check input radius of curv. signs
    if sign(r1) == sign(r2) > 0
        orientation = :left_facing
        convex_radius = r1
        concave_radius = r2
    elseif sign(r1) == sign(r2) < 0
        orientation = :right_facing
        convex_radius = abs(r2)
        concave_radius = abs(r1)
    else
        throw(ArgumentError("Invalid sign combination for r₁ and r₂"))
    end
    check_sag(convex_radius, d)
    check_sag(concave_radius, d)
    # Calculate convex and concave sag
    convex_sag = sag(convex_radius, d)
    concave_sag = sag(concave_radius, d)
    cylinder_l = l - convex_sag + concave_sag
    if cylinder_l ≤ 0
        throw(ErrorException("Lens parameters leads to zero lens edge thickness"))
    end
    # Spawn subshapes
    convex_surface = ConvexSphericalSurfaceSDF(convex_radius, d)
    concave_surface = SphereSDF(concave_radius)
    cylinder = PlanoSurfaceSDF(cylinder_l, d)
    # Move subshapes into position
    if orientation == :left_facing
        translate3d!(cylinder, [0, thickness(convex_surface), 0])
        translate3d!(concave_surface, [0, concave_radius + l, 0])
    elseif orientation == :right_facing
        translate3d!(concave_surface, [0, -concave_radius, 0])
        translate3d!(cylinder, [0, -concave_sag, 0])
        zrotate3d!(convex_surface, π)
        translate3d!(convex_surface, [0, thickness(cylinder) - concave_sag + convex_sag, 0])
    else
        throw(ErrorException("Something went very wrong ... :)"))
    end

    # Return shape
    return MeniscusLensSDF{T, typeof(convex_surface), typeof(concave_surface)}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        convex_surface, cylinder, concave_surface,
        l
    )
end

"""
    MeniscusLensSDF(r1::R1, r2::R2, l::L, d::D, md::MD)

Constructs a positive or negative [`MeniscusLensSDF`](@ref) with:

- `r1`: front surface radis or curvature
- `r2`: back surface radis or curvature
- `l`: lens thickness
- `d`: lens diameter, default value is one inch
- `md`: mechanical lens diameter, must be > d
"""
function MeniscusLensSDF(r1::R1, r2::R2, l::L, d::D, md::MD) where {R1, R2, L, D, MD}
    check_md(d, md)
    # generate ring-less shape
    shape = MeniscusLensSDF(r1, r2, l, d)
    # add an outer ring
    _thickness = thickness(shape.cylinder)
    pos = position(shape.cylinder)
    ring = RingSDF(d / 2, (md - d) / 2, _thickness)
    translate3d!(ring, [0, pos[2] + _thickness / 2, 0])
    shape += ring
    return shape
end

function meniscus_lens_sdf(front_surface::AbstractSurface{T1}, front::AbstractSDF{T2},
        back_surface::AbstractSurface{T3}, back::AbstractSDF{T4},
        center_thickness::Real
) where {T1, T2, T3, T4}
    T = promote_type(T1, T2, T3, T4)

    # Determine orientation: for left‐facing (r₁, r₂ > 0) front is convex and back is concave;
    # for right‐facing (r₁, r₂ < 0) front is concave and back is convex.
    r1, r2 = radius(front_surface), radius(back_surface)
    orientation = if sign(r1) == sign(r2) > 0
        :left_facing
    elseif sign(r1) == sign(r2) < 0
        :right_facing
    else
        throw(ArgumentError("Invalid sign combination for r₁ and r₂"))
    end

    # Compute sag values at the clear aperture:
    # Use d1 for the front and d2 for the back.
    convex_sag = edge_thickness(front_surface, front)
    concave_sag = edge_thickness(back_surface, back)

    cylinder_l = center_thickness - convex_sag + concave_sag
    if cylinder_l ≤ 0
        throw(ErrorException("Lens parameters lead to zero lens edge thickness"))
    end

    # Spawn new sub-shapes.
    if orientation == :left_facing
        front = sdf(front_surface, ForwardLeftMeniscusOrientation())
        back = sdf(back_surface, BackwardLeftMeniscusOrientation())
    else  # right-facing: front is concave, back is convex.
        front = sdf(front_surface, ForwardRightMeniscusOrientation())
        back = sdf(back_surface, BackwardRightMeniscusOrientation())
    end

    # Use effective optical diameter d_mid = min(d1, d2)
    d_mid = min(diameter(front_surface), diameter(back_surface))
    cylinder = PlanoSurfaceSDF(cylinder_l, d_mid)

    # Position the sub-shapes.
    if orientation == :left_facing
        translate3d!(cylinder, [0, thickness(front), 0])
        translate3d!(back, [0, radius(back_surface) + center_thickness, 0])
        convex_shape = front
        concave_shape = back
    else
        translate3d!(back, [0, -abs(radius(front_surface)), 0])
        translate3d!(cylinder, [0, -concave_sag, 0])
        zrotate3d!(front_surface, π)
        translate3d!(front_surface, [0, thickness(cylinder) - concave_sag + convex_sag, 0])
        convex_shape = back
        concave_shape = front
    end

    # Build the composite meniscus lens sdf
    shape = MeniscusLensSDF{T, typeof(convex_shape), typeof(concave_shape)}(
        Matrix{T}(I, 3, 3),
        Matrix{T}(I, 3, 3),
        Point3{T}(0),
        convex_shape,
        cylinder,
        concave_shape,
        center_thickness
    )

    return shape
end
