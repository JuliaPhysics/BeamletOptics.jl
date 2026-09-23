"""
    SquarePlanoMirror2D(edge_length)

Constructs a 2D square plano [`Mirror`](@ref) with a given `edge_length`.
The reflecting surface is normal to the y-axis.

# Inputs

- `edge_length`: the edge length of the square mirror in [m]
"""
function SquarePlanoMirror2D(size::T) where {T <: Real}
    shape = QuadraticFlatMesh(size)
    return Mirror(shape)
end

"""
    RectangularPlanoMirror(width, height, thickness)

Constructs a rectangular plano [`Mirror`](@ref) based on the input dimensions.
The front reflecting surface is normal to the y-axis and lies at the origin.

# Inputs

- `width`:      of the mirror in x-direction [m] 
- `height`:     of the mirror in z-direction [m] 
- `thickness`:  of the mirror in y-direction [m] 
"""
function RectangularPlanoMirror(width::W, height::H, thickness::T) where {W<:Real,H<:Real,T<:Real}
    shape = CuboidMesh(width, thickness, height)
    translate3d!(shape, [
        -width/2,       # x
        0,              # y
        -height/2,      # z
    ])
    set_new_origin3d!(shape)
    return Mirror(shape)
end

"""
    SquarePlanoMirror(width, thickness)

Constructs a square plano [`Mirror`](@ref) with equal width and height.
The front reflecting surface is normal to the y-axis and lies at the origin.
See also [`RectangularPlanoMirror`](@ref).

# Inputs

- `width`: the side length of the square mirror in x- and y-direction [m]
- `thickness`: of the mirror in [m]
"""
function SquarePlanoMirror(width::W, thickness::T) where {W<:Real,T<:Real}
    return RectangularPlanoMirror(width, width, thickness)
end

"""
    RoundPlanoMirror(diameter, thickness; hole_diameter=nothing)

Constructs a round plano [`Mirror`](@ref) with a flat reflecting surface and perfect reflectivity (R = 1).
The reflecting surface is modeled using a [`BeamletOptics.PlanoSurfaceSDF`](@ref).

# Inputs

- `diameter`: mirror diameter \\[m\\]
- `thickness`: mirror substrate thickness \\[m\\]
- `hole_diameter`: diameter of the central through-hole \\[m\\], no hole if `nothing` (default)
"""
function RoundPlanoMirror(diameter::D, thickness::T; hole_diameter::Union{Real, Nothing} = nothing) where {D<:Real,T<:Real}
    shape = PlanoSurfaceSDF(thickness, diameter)
    if hole_diameter === nothing
        return Mirror(shape)
    end
    T_res = typeof(float(thickness))
    return Mirror(_pierce_substrate(shape, diameter, zero(T_res), T_res(thickness), hole_diameter))
end

"""
    RightAnglePrismMirror(leg_length, height)

Constructs a right-angle prism [`Mirror`](@ref) with perfect reflectivity (R = 1).
The reflecting surface is modeled using a [`BeamletOptics.RightAnglePrismSDF`](@ref).

# Inputs

- `leg_length`: edge length in x and y \\[m\\]
- `height`: height in z-axis \\[m\\]
"""
function RightAnglePrismMirror(leg_length::Real, height::Real)
    shape = RightAnglePrismSDF(leg_length, height)
    zrotate3d!(shape, deg2rad(45+180))
    return Mirror(shape)
end
