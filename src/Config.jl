"""
    Config

Global configuration settings for the BeamletOptics package using Preferences.jl.
"""
module Config

using Preferences

# --- Paraxial Invariants ---
const INVARIANT_THRESHOLD = @load_preference("invariant_threshold", 1e-6)
"""
    get_invariant_threshold()

Returns the global configuration threshold for the paraxial optical invariant check.
"""
get_invariant_threshold() = INVARIANT_THRESHOLD

# --- SDF & Ray Marching ---
const SDF_SURFACE_THRESHOLD = @load_preference("sdf_surface_threshold", 1e-9)
const SDF_RAYMARCH_EPS = @load_preference("sdf_raymarch_eps", 1e-13)
const SDF_INSIDE_STEP = @load_preference("sdf_inside_step", 1.0)

"""
    get_sdf_surface_threshold()

Returns the global configuration threshold for the Signed Distance Field (SDF) surface detection.
"""
get_sdf_surface_threshold() = SDF_SURFACE_THRESHOLD

"""
    get_sdf_raymarch_eps()

Returns the global configuration epsilon for the SDF ray marching algorithm.
"""
get_sdf_raymarch_eps() = SDF_RAYMARCH_EPS

"""
    get_sdf_inside_step()

Returns the global configuration step size when inside an SDF volume.
"""
get_sdf_inside_step() = SDF_INSIDE_STEP

# --- Numerical Thresholds ---
const INTERNAL_REFLECTION_THRESHOLD = @load_preference("internal_reflection_threshold",
    1e-6)
const LINE_PLANE_INTERSECTION_THRESHOLD = @load_preference("line_plane_intersection_threshold",
    1e-6)
const ORTHOGONALITY_THRESHOLD = @load_preference("orthogonality_threshold", 1e-10)
const COINCIDENT_BOUNDARY_TOLERANCE = @load_preference("coincident_boundary_tolerance",
    1e-7)
const INDEX_MATCHING_TOLERANCE = @load_preference("index_matching_tolerance", 1e-5)

get_internal_reflection_threshold() = INTERNAL_REFLECTION_THRESHOLD

"""
    get_line_plane_intersection_threshold()

Returns the global configuration threshold for line-plane intersection calculations.
"""
get_line_plane_intersection_threshold() = LINE_PLANE_INTERSECTION_THRESHOLD

"""
    get_orthogonality_threshold()

Returns the global configuration threshold for orthogonality checks (e.g., beam support axes).
"""
get_orthogonality_threshold() = ORTHOGONALITY_THRESHOLD

"""
    get_coincident_boundary_tolerance()

Returns the global configuration threshold for identifying coincident boundaries (e.g., doublet lenses, prism-coating).
"""
get_coincident_boundary_tolerance() = COINCIDENT_BOUNDARY_TOLERANCE

"""
    get_index_matching_tolerance()

Returns the global configuration threshold for comparing refractive indices at interfaces.
"""
get_index_matching_tolerance() = INDEX_MATCHING_TOLERANCE

# --- Tracing Defaults ---
const DEFAULT_R_MAX = @load_preference("default_r_max", 100)
const DEFAULT_DEPTH_MAX = @load_preference("default_depth_max", typemax(Int))
const DEFAULT_POWER_CUTOFF = @load_preference("default_power_cutoff", 0.0)

"""
    get_default_r_max()

Returns the default maximum number of segments for ray tracing.
"""
get_default_r_max() = DEFAULT_R_MAX

"""
    get_default_depth_max()

Returns the default maximum depth for recursive ray tracing (e.g., reflections/refractions).
"""
get_default_depth_max() = DEFAULT_DEPTH_MAX

"""
    get_default_power_cutoff()

Returns the default power cutoff threshold below which sub-beams/split paths are dropped to prevent infinite branching.
"""
get_default_power_cutoff() = DEFAULT_POWER_CUTOFF

# --- Physical Defaults ---
const DEFAULT_WAVELENGTH = @load_preference("default_wavelength", 1e-6)
const DEFAULT_WAIST = @load_preference("default_waist", 1e-3)
const DEFAULT_POWER = @load_preference("default_power", 1e-3)

"""
    get_default_wavelength()

Returns the default wavelength in meters.
"""
get_default_wavelength() = DEFAULT_WAVELENGTH

"""
    get_default_waist()

Returns the default beam waist radius in meters.
"""
get_default_waist() = DEFAULT_WAIST

"""
    get_default_power()

Returns the default beam total power in Watts.
"""
get_default_power() = DEFAULT_POWER

"""
    set_preference!(key::String, val; persistent=true)

Helper to set preferences. Default is persistent.
"""
function set_preference!(key::String, val)
    @set_preferences!(key=>val)
    @info "Preference '$key' updated to $val. Please restart Julia for this to take effect."
end

# Keep the old setter for compatibility
"""
    set_invariant_threshold!(val::Real)

Sets the global configuration threshold for the paraxial optical invariant check.
This preference is persistent. Please restart Julia for this to take effect.
"""
set_invariant_threshold!(val::Real) = set_preference!("invariant_threshold", Float64(val))

end # module
