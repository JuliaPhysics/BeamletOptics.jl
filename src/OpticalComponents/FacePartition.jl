# Fallback for transposed orientation
transposed_orientation(shape::AbstractShape) = transpose(orientation(shape))

"""
    world_to_local(shape::AbstractShape, point::AbstractVector)

Transforms coordinate `point` from world space to the local coordinate system of `shape`.
"""
function world_to_local(shape::AbstractSDF, point::AbstractVector)
    return _world_to_sdf(shape, point)
end

function world_to_local(shape::AbstractShape, point::AbstractVector)
    T = transposed_orientation(shape)
    return T * (point - position(shape))
end

"""
    face_id(shape::AbstractShape, local_n::AbstractVector)::Symbol

Partition the surface normal space of a shape into named faces (e.g., `:front`, `:back`, `:side`, `:hypotenuse`).
Returns a `Symbol` corresponding to the face that contains the given `local_n` vector (in local shape coordinates).
"""
function face_id(shape::AbstractShape, local_n::AbstractVector)
    # Default fallback for rotationally symmetric/general shapes:
    # y-axis is the optical axis in BMO... *sigh* ... Why?
    # Front is pointing in -y direction.
    # Back is pointing in +y direction.
    ny = local_n[2]
    if ny < -0.1
        return :front
    elseif ny > 0.1
        return :back
    else
        return :side
    end
end

# Specialization for RightAnglePrismSDF
function face_id(shape::RightAnglePrismSDF, local_n::AbstractVector)
    # Face normals:
    # - Hypotenuse: [1.0, 1.0, 0.0] / √2
    # - Leg 1: [-1.0, 0.0, 0.0]
    # - Leg 2: [0.0, -1.0, 0.0]
    # - Zpos: [0.0, 0.0, 1.0]
    # - Zneg: [0.0, 0.0, -1.0]
    d_hyp = dot(local_n, SVector(1.0, 1.0, 0.0) / √2)
    d_leg1 = dot(local_n, SVector(-1.0, 0.0, 0.0))
    d_leg2 = dot(local_n, SVector(0.0, -1.0, 0.0))
    d_zpos = dot(local_n, SVector(0.0, 0.0, 1.0))
    d_zneg = dot(local_n, SVector(0.0, 0.0, -1.0))

    best_face = argmax((d_hyp, d_leg1, d_leg2, d_zpos, d_zneg))
    return (:hypotenuse, :leg1, :leg2, :zpos, :zneg)[best_face]
end
