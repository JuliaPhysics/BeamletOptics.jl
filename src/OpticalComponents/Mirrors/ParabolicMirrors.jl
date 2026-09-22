"""
    OffAxisParabolicMirror(rfl, diameter; angle=90, thickness=nothing, hole_diameter=nothing, hole_axis=:collimated)

Constructs an Off-Axis Parabolic (OAP) [`Mirror`](@ref) from:

# Inputs

- `rfl`:            Reflected Focal Length (distance from aperture center to focus) [m]
- `diameter`:       Mirror aperture diameter [m]
- `angle`:          Deflection angle in degrees (default: 90°)
- `thickness`:      Substrate thickness [m], calculated automatically to ensure solid backing if `nothing` (default)
- `hole_diameter`:  Diameter of the through-hole [m], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
- `hole_axis`:      Orientation of the through-hole. Options:
                    - `:collimated` (default): parallel to the collimated beam (local y-axis / substrate normal), centered at the aperture center `(0, 0, 0)`.
                    - `:focused`: angled towards the parent paraboloid focus `(-x_off, x_off^2/(4f) - f, 0)` in the segment frame, equivalently `(-rfl*sind(angle), -rfl*cosd(angle), 0)` with `angle` in degrees, passing through the aperture center `(0, 0, 0)` (e.g. for collinear pump-probe beams).
"""
function OffAxisParabolicMirror(
        rfl::Real,
        diameter::Real;
        angle::Real = 90,
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing,
        hole_axis::Symbol = :collimated
    )
    T = float(promote_type(typeof(rfl), typeof(diameter), typeof(angle), typeof(thickness === nothing ? 0.0 : thickness)))
    angle_rad = deg2rad(angle)

    f = T(rfl * (cos(angle_rad / 2)^2))
    x_off = T(rfl * sin(angle_rad))

    r_max = T(diameter / 2)
    sag_max = abs(-(((r_max + x_off)^2 - x_off^2) / (4 * f)))

    t = thickness === nothing ? max(T(diameter / 2), sag_max + T(10e-3)) : T(thickness)

    oap_sdf = ConicSDF(2f, -one(T), x_off, T(diameter), t)
    if hole_diameter === nothing
        return Mirror(oap_sdf)
    end

    hd = T(hole_diameter)
    _check_hole_diameter(hd, T(diameter))

    if hole_axis === :collimated
        return Mirror(_pierce_substrate(oap_sdf, diameter, -sag_max, t, hole_diameter))
    elseif hole_axis === :focused
        focus_vec = Point3(-x_off, -T(rfl * cos(angle_rad)), zero(T))
        u_dir = normalize(focus_vec)
        margin = T(10e-3)
        span = max(T(diameter), t + sag_max) + 2margin
        bore = CylinderSDF(hd / 2, span)
        align3d!(bore, u_dir)
        return Mirror(oap_sdf - bore)
    else
        throw(ArgumentError("hole_axis must be :collimated or :focused, got :$hole_axis"))
    end
end

"""
    ParabolicMirror(f, diameter; thickness=nothing, hole_diameter=nothing)

Constructs an on-axis parabolic [`Mirror`](@ref) with focal length `f`.
The vertex of the concave reflecting surface lies at the origin, the mirror opens towards the negative y-axis
and its focus lies at `(0, -f, 0)`. The shape is an [`OffAxisParaboloidSDF`](@ref) without off-axis offset.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (a Cassegrain primary).

!!! warning "Reflective bore wall"
    Rays that graze into the hole will reflect off its wall rather than being absorbed.

# Inputs

- `f`:              Focal length \\[m\\]
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], rim sag + 10 mm if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function ParabolicMirror(
        f::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    T = float(promote_type(typeof(f), typeof(diameter), typeof(thickness === nothing ? 0.0 : thickness)))
    sag_max = T(diameter / 2)^2 / (4 * T(f))
    t = thickness === nothing ? sag_max + T(10e-3) : T(thickness)
    substrate = ConicSDF(2 * T(f), -one(T), zero(T), T(diameter), t)
    if hole_diameter === nothing
        return Mirror(substrate)
    end
    return Mirror(_pierce_substrate(substrate, diameter, -sag_max, t, hole_diameter))
end
