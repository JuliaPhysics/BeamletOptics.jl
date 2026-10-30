"""
    OffAxisParabolicMirror(rfl, diameter; angle=90, thickness=nothing, hole_diameter=nothing, hole_axis=:collimated)

Constructs an Off-Axis Parabolic (OAP) [`Mirror`](@ref) from:

# Arguments

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

    # Parent focus in the segment frame
    focus = Point3(-x_off, -T(rfl * cos(angle_rad)), zero(T))
    return build_oap(hole_axis, f, x_off, diameter, t, sag_max, focus, hole_diameter)
end

function build_oap(::Symbol, f::T, x_off::T, diameter::Real, t::T, ::T, ::Point3{T},
        ::Nothing) where {T <: Real}
    return OffAxisConicMirror(2f, -1, x_off, diameter; thickness = t)
end

function build_oap(hole_axis::Symbol, f::T, x_off::T, diameter::Real, t::T, sag_max::T,
        focus::Point3{T}, hole_diameter::Real) where {T <: Real}
    if hole_axis === :collimated
        return OffAxisConicMirror(2f, -1, x_off, diameter; thickness = t, hole_diameter)
    elseif hole_axis === :focused
        oap_sdf = ConicSDF(2f, -one(T), x_off, T(diameter), t)
        hd = T(hole_diameter)
        check_hole_diameter(hd, T(diameter))
        u_dir = normalize(focus)
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

Constructs an on-axis parabolic [`Mirror`](@ref) with focal length `f`, i.e. a
[`ConicMirror`](@ref) with `R = 2f` and `k = -1`. The vertex of the reflecting surface lies at
the origin and the focus at `(0, -f, 0)`. For `f > 0` the mirror is concave and opens towards
the negative y-axis; for `f < 0` it is convex and its focus is virtual, behind the mirror. The
shape is an [`OffAxisParaboloidSDF`](@ref) without off-axis offset.

If `hole_diameter` is given, a cylindrical bore centred on the optical axis (the local
+y-axis) is subtracted from the substrate, piercing it completely (a Cassegrain primary).

!!! warning "Reflective bore wall"
    Rays that graze into the hole will reflect off its wall rather than being absorbed.

# Arguments

- `f`:              Focal length \\[m\\]; `f > 0` concave, `f < 0` convex
- `diameter`:       Mirror aperture diameter \\[m\\]
- `thickness`:      Substrate thickness \\[m\\], `|rim sag| + 10 mm` if `nothing` (default)
- `hole_diameter`:  Diameter of the central through-hole \\[m\\], no hole if `nothing` (default). Must satisfy `0 < hole_diameter < diameter`.
"""
function ParabolicMirror(
        f::Real,
        diameter::Real;
        thickness::Union{Real, Nothing} = nothing,
        hole_diameter::Union{Real, Nothing} = nothing
    )
    return ConicMirror(2f, -1, diameter; thickness, hole_diameter)
end
