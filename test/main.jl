using SCDI
using LinearAlgebra
using FileIO, GLMakie

mutable struct Geometry
    vertices::Matrix
    faces::Matrix
    pos::Vector
    dir::Vector
    function Geometry(data)
        # Read data from shitty mesh format to matrix format
        numEl = length(data)
        vertices = Matrix{Float64}(undef, numEl*3, 3)
        faces = Matrix{Int}(undef, numEl, 3)
        for i = 1:numEl
            for j = 1:3
                vertices[3(i-1)+j, :] = data[i][j]
                faces[i,j] = 3(i-1)+j
            end
        end
        # Initialize geometry at origin with orientation [1,0,0] scaled to mm
        new(vertices*1e-3, faces, [0,0,0], [1,0,0])
    end
end

function translate3d!(object::Geometry, offset::Vector)
    object.pos += offset
    object.vertices .+= offset'
end

function xrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = SCDI.rotate3d([1,0,0], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

function yrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = SCDI.rotate3d([0,1,0], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

function zrotate3d!(object::Geometry, θ)
    # Calculate rotation matrix
    R = SCDI.rotate3d([0,0,1], θ)
    # Translate geometry to origin, rotate, retranslate
    object.vertices = (object.vertices .- object.pos')*R .+ object.pos'
end

# face i can be accessed via matrix as cube.vertices[cube.faces[i,:],:]
# face coordinates must be global!
# MT algorithm is fast due to breakout conditions
function intersect3d(face::Matrix, ray::SCDI.Ray; kϵ=1e-5)
    V1 = face[1,:]
    V2 = face[2,:]
    V3 = face[3,:]

    E1 = V2-V1
    E2 = V3-V1
    Pv = cross(ray.dir, E2)
    Det = dot(E1, Pv)
    invDet = 1/Det
    # Check if ray is backfacing or missing face
    if (Det < kϵ) & (abs(Det) < kϵ)
        return Inf
    end
    # Compute normalized u and reject if less than 0 or greater than 1
    Tv = ray.pos-V1
    u = dot(dot(Tv, Pv), invDet)
    if (u < 0) | (u > 1)
        return Inf
    end
    # Compute normalized v and reject if less than 0 or greater than 1
    Qv = cross(Tv, E1)
    v = dot(dot(ray.dir, Qv), invDet)
    if (v < 0) | (v > 1)
        return Inf
    end
    # Compute t and return positive intersection (type def. for t to avoid Any)
    t::Float64 = dot(dot(E2, Qv), invDet)
    return t
end

function intersect3d(object::Geometry, ray::SCDI.Ray)
    numEl = size(object.faces)[1]
    t0 = Inf
    for i = 1:numEl
        face = object.vertices[object.faces[i,:],:]
        t = intersect3d(face, ray)
        if t < t0
            t0 = t
        end
    end
    return t0
end

cube = Geometry(load("test\\cube.stl"))
translate3d!(cube, [0.01,-0.01,-0.01])
zrotate3d!(cube, π/4)
xrotate3d!(cube, π/4)

f = Figure()
ax = f[1,1] = Axis3(f,aspect = (1,1,1))
xlims!(ax, (-0.025, 0.025))
ylims!(ax, (-0.025, 0.025))
zlims!(ax, (-0.025, 0.025))
mesh!(cube.vertices, cube.faces)
f

ray = SCDI.Ray([-1,0,0], [1,0,0])
t = intersect3d(cube, ray)

lines!([-1,-1+t],[0,0],[0,0], color=:red)
f
