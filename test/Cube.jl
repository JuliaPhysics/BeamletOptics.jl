### Helper functions for testing

function Cube(scale::Real; T = Float64)
    vertices = [0 0 0
        1 0 0
        1 1 0
        0 1 0
        0 1 1
        1 1 1
        1 0 1
        0 0 1]
    faces = [1 3 2
        1 4 3
        3 4 5
        3 5 6
        2 3 6
        2 6 7
        1 8 5
        1 5 4
        6 5 8
        6 8 7
        1 7 8
        1 2 7]
    return SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

struct ReflectiveCube{S <: SCDI.AbstractShape} <: SCDI.AbstractReflectiveOptic
    id::UUID
    shape::S
end

function RetroMesh(scale::Real; T = Float64)
    vertices = [0 0 0
        1 0 0
        0 1 0
        0 0 1]
    faces = [1 3 2
        1 4 3
        1 2 4]
    return SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
end

struct RetroReflector{S <: SCDI.AbstractShape} <: SCDI.AbstractReflectiveOptic
    id::UUID
    shape::S
end

RetroReflector(scale) = RetroReflector(uuid4(), RetroMesh(scale))

function PlanoMirror(scale::T) where {T <: Real}
    sz = 0.5
    vertices = [sz 0 sz
        sz 0 -sz
        -sz 0 -sz
        -sz 0 sz]
    faces = [1 2 4
        2 3 4]
    shape = SCDI.Mesh{T}(uuid4(),
        vertices .* scale,
        faces,
        Matrix{T}(I, 3, 3),
        T.([0, 0, 0]),
        scale)
    return SCDI.Mirror(uuid4(), shape)
end

using StaticArrays

struct ObjectGroup{T} <: SCDI.AbstractObjectGroup
    dir::Matrix{T}
    center::Vector{T}
    objects::Vector{SCDI.AbstractObject}
end

SCDI.position(group::ObjectGroup) = group.center
SCDI.position!(group::ObjectGroup, pos) = (group.center .= pos)

SCDI.orientation(group::ObjectGroup) = group.dir
SCDI.orientation!(group::ObjectGroup, dir) = (group.dir .= dir)

function ObjectGroup(v::AbstractArray{<:SCDI.AbstractObject}, T = Float64)
    ObjectGroup{T}(Matrix{T}(I, 3, 3), T.([0, 0, 0]), v)
end

SCDI.objects(group::ObjectGroup) = group.objects

function SCDI.translate3d!(group::ObjectGroup, offset)
    # Translate tracking vector
    SCDI.position!(group, SCDI.position(group) .+ offset)
    # Recursively translate all subgroups
    for object in SCDI.objects(group)
        SCDI.translate3d!(object, offset)
    end
    return nothing
end

function SCDI.rotate3d!(group::ObjectGroup, axis, θ)
    R = SCDI.rotate3d(axis, θ)
    # Update group orientation
    SCDI.orientation!(group, SCDI.orientation(group) * R)
    # Recursively rotate all subgroups and objects
    for object in SCDI.objects(group)
        SCDI.rotate3d!(object, axis, θ)
        v = SCDI.position(object) - SCDI.position(group)
        # Translate group around pivot point
        v = (R * v) - v
        SCDI.translate3d!(object, v)
    end
    return nothing
end

function SCDI.xrotate3d!(group::ObjectGroup{T}, θ) where {T}
    SCDI.rotate3d!(group, @SVector(T[one(T), zero(T), zero(T)]), θ)
end
function SCDI.yrotate3d!(group::ObjectGroup{T}, θ) where {T}
    SCDI.rotate3d!(group, @SVector(T[zero(T), one(T), zero(T)]), θ)
end
function SCDI.zrotate3d!(group::ObjectGroup{T}, θ) where {T}
    SCDI.rotate3d!(group, @SVector(T[zero(T), zero(T), one(T)]), θ)
end

Base.show(::IO, ::MIME"text/plain", group::ObjectGroup) = print_tree(group)
