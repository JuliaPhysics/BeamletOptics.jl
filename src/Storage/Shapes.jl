#=
Storage of shapes, see `save_setup`.

Leaf SDFs are stored by their constructor parameters plus their world `pos` and `dir`. To load
them, the leaf is constructed at the origin and then moved with `orientation!` and `position!`,
which keep the cached `transposed_dir` in sync. Derived fields (e.g. `sag`, `w`, `max_sag`) are
recomputed by the constructors.

A `UnionSDF` stores its own pose and its sub-SDFs at their world poses. A `MeniscusLensSDF`
stores its own pose and its sub-SDFs at their poses relative to the lens. A `Mesh` stores its
vertices (already at world pose) and faces as arrays, which become binary assets if large.
=#

function _check_shape_eltype(s::AbstractShape{T}) where {T}
    T === Float64 || throw(ArgumentError("only Float64 shapes can be stored, got $(typeof(s))"))
    return nothing
end

_float_vector(x) = Float64[Float64(v) for v in x]

# Pose and parameters of a leaf SDF
function _sdf_storage(s::AbstractSDF, params::Pair...)
    _check_shape_eltype(s)
    d = Dict{String, Any}("pos" => encode_vec3(position(s)), "dir" => encode_mat3(orientation(s)))
    for (k, v) in params
        d[k] = v
    end
    return d
end

# Moves a leaf SDF constructed at the origin to its stored pose
function _set_pose!(s::AbstractSDF, d)
    orientation!(s, decode_mat3(d["dir"]))
    position!(s, decode_vec3(d["pos"]))
    return s
end

#=
Primitives
=#

# The constructor expects the full edge lengths, the struct holds half of them
to_storage(s::BoxSDF, ctx) = _sdf_storage(s, "size" => encode_vec3(2 .* s.dimensions))
function from_storage(::Type{BoxSDF}, d, ctx)
    x, y, z = decode_vec3(d["size"])
    return _set_pose!(BoxSDF(x, y, z), d)
end

to_storage(s::CylinderSDF, ctx) = _sdf_storage(s, "radius" => s.radius, "height" => s.height)
from_storage(::Type{CylinderSDF}, d, ctx) =
    _set_pose!(CylinderSDF(Float64(d["radius"]), Float64(d["height"])), d)

to_storage(s::CutSphereSDF, ctx) = _sdf_storage(s, "radius" => s.radius, "height" => s.height)
from_storage(::Type{CutSphereSDF}, d, ctx) =
    _set_pose!(CutSphereSDF(Float64(d["radius"]), Float64(d["height"])), d)

# The constructor expects inner radius, width and thickness, the struct holds the ring center
# radius and half of width and thickness
to_storage(s::RingSDF, ctx) = _sdf_storage(s, "inner_radius" => s.inner_radius - s.hwidth,
    "width" => 2 * s.hwidth, "thickness" => 2 * s.hthickness)
from_storage(::Type{RingSDF}, d, ctx) = _set_pose!(
    RingSDF(Float64(d["inner_radius"]), Float64(d["width"]), Float64(d["thickness"])), d)

to_storage(s::RightAnglePrismSDF, ctx) =
    _sdf_storage(s, "leg_length" => 2 * s.dimensions[1], "height" => 2 * s.dimensions[3])
from_storage(::Type{RightAnglePrismSDF}, d, ctx) =
    _set_pose!(RightAnglePrismSDF(Float64(d["leg_length"]), Float64(d["height"])), d)

#=
Spherical lens surfaces
=#

to_storage(s::PlanoSurfaceSDF, ctx) =
    _sdf_storage(s, "thickness" => s.thickness, "diameter" => s.diameter)
from_storage(::Type{PlanoSurfaceSDF}, d, ctx) =
    _set_pose!(PlanoSurfaceSDF(Float64(d["thickness"]), Float64(d["diameter"])), d)

# The orientation of a sphere is fixed
function to_storage(s::SphereSDF, ctx)
    _check_shape_eltype(s)
    return Dict{String, Any}("pos" => encode_vec3(position(s)), "radius" => s.radius)
end
function from_storage(::Type{SphereSDF}, d, ctx)
    s = SphereSDF(Float64(d["radius"]))
    position!(s, decode_vec3(d["pos"]))
    return s
end

to_storage(s::ConcaveSphericalSurfaceSDF, ctx) =
    _sdf_storage(s, "radius" => s.radius, "diameter" => s.diameter)
from_storage(::Type{ConcaveSphericalSurfaceSDF}, d, ctx) =
    _set_pose!(ConcaveSphericalSurfaceSDF(Float64(d["radius"]), Float64(d["diameter"])), d)

to_storage(s::ConvexSphericalSurfaceSDF, ctx) =
    _sdf_storage(s, "radius" => s.radius, "diameter" => s.diameter)
from_storage(::Type{ConvexSphericalSurfaceSDF}, d, ctx) =
    _set_pose!(ConvexSphericalSurfaceSDF(Float64(d["radius"]), Float64(d["diameter"])), d)

# The sub-SDFs are evaluated in the lens frame, i.e. their poses are relative to the lens
to_storage(s::MeniscusLensSDF, ctx) = _sdf_storage(s,
    "convex" => encode(s.convex, ctx), "cylinder" => encode(s.cylinder, ctx),
    "concave" => encode(s.concave, ctx), "thickness" => s.thickness)
function from_storage(::Type{MeniscusLensSDF}, d, ctx)
    convex = decode(d["convex"], ctx)
    cylinder = decode(d["cylinder"], ctx)
    concave = decode(d["concave"], ctx)
    cylinder isa PlanoSurfaceSDF{Float64} ||
        throw(ArgumentError("the cylinder of a MeniscusLensSDF must be a PlanoSurfaceSDF, got $(typeof(cylinder))"))
    I3 = SMatrix{3, 3, Float64, 9}(I)
    s = MeniscusLensSDF{Float64, typeof(convex), typeof(concave)}(
        I3, I3, Point3{Float64}(0), convex, cylinder, concave, Float64(d["thickness"]))
    return _set_pose!(s, d)
end

#=
Aspheric and cylindric lens surfaces
=#

to_storage(s::ConvexAsphericalSurfaceSDF, ctx) = _sdf_storage(s,
    "coefficients" => collect(Float64, s.coefficients), "radius" => s.radius,
    "conic_constant" => s.conic_constant, "diameter" => s.diameter)
from_storage(::Type{ConvexAsphericalSurfaceSDF}, d, ctx) = _set_pose!(
    ConvexAsphericalSurfaceSDF(_float_vector(d["coefficients"]), Float64(d["radius"]),
        Float64(d["conic_constant"]), Float64(d["diameter"])), d)

to_storage(s::ConcaveAsphericalSurfaceSDF, ctx) = _sdf_storage(s,
    "coefficients" => collect(Float64, s.coefficients), "radius" => s.radius,
    "conic_constant" => s.conic_constant, "diameter" => s.diameter,
    "mechanical_diameter" => s.mechanical_diameter)
from_storage(::Type{ConcaveAsphericalSurfaceSDF}, d, ctx) = _set_pose!(
    ConcaveAsphericalSurfaceSDF(_float_vector(d["coefficients"]), Float64(d["radius"]),
        Float64(d["conic_constant"]), Float64(d["diameter"]), Float64(d["mechanical_diameter"])), d)

to_storage(s::ConvexCylinderSDF, ctx) =
    _sdf_storage(s, "radius" => s.radius, "diameter" => s.diameter, "height" => s.height)
from_storage(::Type{ConvexCylinderSDF}, d, ctx) = _set_pose!(
    ConvexCylinderSDF(Float64(d["radius"]), Float64(d["diameter"]), Float64(d["height"])), d)

to_storage(s::ConcaveCylinderSDF, ctx) =
    _sdf_storage(s, "radius" => s.radius, "diameter" => s.diameter, "height" => s.height)
from_storage(::Type{ConcaveCylinderSDF}, d, ctx) = _set_pose!(
    ConcaveCylinderSDF(Float64(d["radius"]), Float64(d["diameter"]), Float64(d["height"])), d)

_acylinder_storage(s, ctx) = _sdf_storage(s, "radius" => s.radius, "diameter" => s.diameter,
    "height" => s.height, "conic_constant" => s.conic_constant,
    "coefficients" => collect(Float64, s.coefficients))
_acylinder_args(d) = (Float64(d["radius"]), Float64(d["diameter"]), Float64(d["height"]),
    Float64(d["conic_constant"]), _float_vector(d["coefficients"]))

to_storage(s::AconvexCylinderSDF, ctx) = _acylinder_storage(s, ctx)
from_storage(::Type{AconvexCylinderSDF}, d, ctx) = _set_pose!(AconvexCylinderSDF(_acylinder_args(d)...), d)

to_storage(s::AconcaveCylinderSDF, ctx) = _acylinder_storage(s, ctx)
from_storage(::Type{AconcaveCylinderSDF}, d, ctx) = _set_pose!(AconcaveCylinderSDF(_acylinder_args(d)...), d)

to_storage(s::OffAxisParaboloidSDF, ctx) = _sdf_storage(s,
    "f" => s.f, "x_off" => s.x_off, "diameter" => s.diameter, "thickness" => s.thickness)
from_storage(::Type{OffAxisParaboloidSDF}, d, ctx) = _set_pose!(
    OffAxisParaboloidSDF(Float64(d["f"]), Float64(d["x_off"]), Float64(d["diameter"]),
        Float64(d["thickness"])), d)

#=
Unions
=#

# The sub-SDFs are stored at their world poses. The pose of the union is set without moving them,
# i.e. not with `translate3d!` or `rotate3d!`.
to_storage(s::UnionSDF, ctx) = _sdf_storage(s, "sdfs" => [encode(sub, ctx) for sub in s.sdfs])
function from_storage(::Type{UnionSDF}, d, ctx)
    sdfs = [decode(sub, ctx) for sub in d["sdfs"]]
    for sub in sdfs
        sub isa AbstractSDF{Float64} ||
            throw(ArgumentError("UnionSDF members must be Float64 SDFs, got $(typeof(sub))"))
    end
    return _set_pose!(UnionSDF{Float64}(sdfs...), d)
end

#=
Meshes
=#

# The vertices are stored at their world pose, `pos`, `dir` and `scale` only track it
function to_storage(m::Mesh, ctx)
    _check_shape_eltype(m)
    return Dict{String, Any}(
        "vertices" => encode_array(m.vertices, ctx; stem = "mesh-vertices"),
        "faces" => encode_array(collect(Int64, m.faces), ctx; stem = "mesh-faces"),
        "pos" => encode_vec3(m.pos), "dir" => encode_mat3(m.dir), "scale" => m.scale)
end
function from_storage(::Type{Mesh}, d, ctx)
    vertices = decode_array(d["vertices"], ctx)
    faces = decode_array(d["faces"], ctx)
    (vertices isa Matrix{Float64} && size(vertices, 2) == 3) ||
        throw(ArgumentError("mesh vertices must be a Float64 matrix with 3 columns"))
    (faces isa Matrix{Int64} && size(faces, 2) == 3) ||
        throw(ArgumentError("mesh faces must be an Int64 matrix with 3 columns"))
    return Mesh{Float64}(vertices, faces, decode_mat3(d["dir"]), decode_vec3(d["pos"]),
        Float64(d["scale"]))
end

for T in (BoxSDF, CylinderSDF, CutSphereSDF, RingSDF, RightAnglePrismSDF, PlanoSurfaceSDF,
        SphereSDF, ConcaveSphericalSurfaceSDF, ConvexSphericalSurfaceSDF, MeniscusLensSDF,
        ConvexAsphericalSurfaceSDF, ConcaveAsphericalSurfaceSDF, ConvexCylinderSDF,
        ConcaveCylinderSDF, AconvexCylinderSDF, AconcaveCylinderSDF, OffAxisParaboloidSDF,
        UnionSDF, Mesh)
    register_storage_type!(T, string(nameof(T)))
end
