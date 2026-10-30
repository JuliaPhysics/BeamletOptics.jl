"""
    Config

Global configuration settings for the BeamletOptics package using Preferences.jl.
"""
module Config

using Preferences

# --- Paraxial Invariants ---
const invariant_threshold = @load_preference("invariant_threshold", 1e-6)
"""
    get_invariant_threshold()

Returns the global configuration threshold for the paraxial optical invariant check.
"""
get_invariant_threshold() = invariant_threshold


# --- SDF & Ray Marching ---
const sdf_surface_threshold = @load_preference("sdf_surface_threshold", 1e-9)
const sdf_raymarch_eps = @load_preference("sdf_raymarch_eps", 1e-10)
const sdf_miss_distance = @load_preference("sdf_miss_distance", 1e15)

"""
    get_sdf_surface_threshold()

Returns the global configuration threshold for the Signed Distance Field (SDF) surface detection.
"""
get_sdf_surface_threshold() = sdf_surface_threshold

"""
    get_sdf_raymarch_eps()

Returns the global configuration epsilon for the SDF ray marching algorithm.
"""
get_sdf_raymarch_eps() = sdf_raymarch_eps

"""
    get_sdf_miss_distance()

Returns the path length in \\[m\\] after which a ray that is marched towards an SDF shape counts as
a miss: a shape farther from the start of the ray than this is not hit. Defaults to `1e15`
(configurable via the `sdf_miss_distance` preference).
"""
get_sdf_miss_distance() = sdf_miss_distance


# --- Coincident boundaries ---
const coincidence_tolerance = @load_preference("coincidence_tolerance", 1e-9)

"""
    get_coincidence_tolerance()

Returns the path length in [m] within which the hits of a ray with several surfaces count as one
point, e.g. the faces of two cemented lenses or a coating on a substrate. Gaps below this length are
taken as contact. The default of 1 nm lies above the accuracy of the intersection algorithms and far
below any gap that matters optically.
"""
get_coincidence_tolerance() = coincidence_tolerance


# --- Numerical Thresholds ---
const internal_reflection_threshold = @load_preference("internal_reflection_threshold", 1e-6)
const line_plane_intersection_threshold = @load_preference("line_plane_intersection_threshold", 1e-6)
const orthogonality_threshold = @load_preference("orthogonality_threshold", 1e-10)

"""
    get_internal_reflection_threshold()

Returns the global configuration threshold for total internal reflection detection.
"""
get_internal_reflection_threshold() = internal_reflection_threshold

"""
    get_line_plane_intersection_threshold()

Returns the global configuration threshold for line-plane intersection calculations.
"""
get_line_plane_intersection_threshold() = line_plane_intersection_threshold

"""
    get_orthogonality_threshold()

Returns the global configuration threshold for orthogonality checks (e.g., beam support axes).
"""
get_orthogonality_threshold() = orthogonality_threshold


# --- Tracing Defaults ---
const default_r_max = @load_preference("default_r_max", 100)
const default_depth_max = @load_preference("default_depth_max", 100)

"""
    get_default_r_max()

Returns the default maximum number of segments for ray tracing.
"""
get_default_r_max() = default_r_max

"""
    get_default_depth_max()

Returns the default maximum depth for recursive ray tracing (e.g., reflections/refractions).
Defaults to 100 (configurable via the `default_depth_max` preference).
"""
get_default_depth_max() = default_depth_max


# --- Physical Defaults ---
const default_wavelength = @load_preference("default_wavelength", 1e-6)
const default_waist = @load_preference("default_waist", 1e-3)
const default_power = @load_preference("default_power", 1e-3)

"""
    get_default_wavelength()

Returns the default wavelength in meters.
"""
get_default_wavelength() = default_wavelength

"""
    get_default_waist()

Returns the default beam waist radius in meters.
"""
get_default_waist() = default_waist

"""
    get_default_power()

Returns the default beam total power in Watts.
"""
get_default_power() = default_power


# --- Progress Bars ---
const progress_threshold = @load_preference("progress_threshold", 5.0)

"""
    get_progress_threshold()

Returns the runtime in seconds after which long-running calls (`solve_system!` on a beam group,
`electric_field`, `intensity` and `optical_power` on a detector) show a progress bar in the terminal.
Defaults to 5 s (configurable via the `progress_threshold` preference, see
[`set_progress_threshold!`](@ref)). `Inf` disables all progress bars.
"""
get_progress_threshold() = progress_threshold


"""
    set_preference!(key::String, val; persistent=true)

Helper to set preferences. Default is persistent.
"""
function set_preference!(key::String, val)
    @set_preferences!(key => val)
    @info "Preference '$key' updated to $val. Please restart Julia for this to take effect."
end

# Keep the old setter for compatibility
"""
    set_invariant_threshold!(val::Real)

Sets the global configuration threshold for the paraxial optical invariant check.
This preference is persistent. Please restart Julia for this to take effect.
"""
set_invariant_threshold!(val::Real) = set_preference!("invariant_threshold", Float64(val))

"""
    set_progress_threshold!(val::Real)

Sets the runtime in seconds after which long-running calls show a progress bar, see
[`get_progress_threshold`](@ref). Use `Inf` to disable progress bars globally.
This preference is persistent. Please restart Julia for this to take effect.
"""
set_progress_threshold!(val::Real) = set_preference!("progress_threshold", Float64(val))


end # module
