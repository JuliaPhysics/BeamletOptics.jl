using BeamletOptics
using StaticArrays
using LinearAlgebra

# ==============================================================================
# Coatings Demo
# ==============================================================================
# This script demonstrates the usage of Optical Coatings in BeamletOptics.jl.
# We cover:
# 1. Simple Coatings (AR / Mirror)
# 2. Multilayer Coatings (TMM)
# 3. Spatially Varying Coatings (Front vs Back)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Simple Coated Lens
# ------------------------------------------------------------------------------
println("--- 1. Simple AR Coated Lens ---")
# Create a standard singlet lens with wide-band AR coating (0.5% R)
ar_coating = SimpleCoating(0.005) # 0.5% Reflection, 99.5% Transmission

lens_simple = Lens(
    SphericalSurface(50e-3, 25e-3),
    SphericalSurface(-50e-3, 25e-3),
    5e-3,
    x -> 1.5,   # N-BK7 appox
    ar_coating  # Apply to all surfaces
)
println("Created lens with simple AR coating: $(lens_simple.coating)")

# ------------------------------------------------------------------------------
# 2. Multilayer Coating (Quarter-Wave AR)
# ------------------------------------------------------------------------------
println("\n--- 2. Multilayer AR Coating (TMM) ---")
# Design a single layer MgF2 AR coating centered at 550nm
λ_design = 550e-9
n_MgF2 = 1.38
d_layer = λ_design / (4 * n_MgF2)

layer = ThinFilmLayer(n_MgF2, d_layer)
multilayer_ar = MultilayerCoating([layer])

# Check reflectivity at design wavelength
n_air = 1.0
n_glass = 1.5
rs, _, _, _ = BeamletOptics.coating_coefficients(
    multilayer_ar, n_air, n_glass, λ_design, 0.0)
R_design = abs2(rs)
println("Reflectivity at $(λ_design*1e9)nm: $(round(R_design*100, digits=4))%")

# Check reflectivity off-design
λ_off = 800e-9
rs_off, _, _, _ = BeamletOptics.coating_coefficients(
    multilayer_ar, n_air, n_glass, λ_off, 0.0)
println("Reflectivity at $(λ_off*1e9)nm:  $(round(abs2(rs_off)*100, digits=4))%")

# ------------------------------------------------------------------------------
# 3. Spatially Varying Coating (Front AR, Back HR)
# ------------------------------------------------------------------------------
println("\n--- 3. Spatially Varying Coating ---")
# Create a Lens that acts as a partial reflector:
# Front Surface: AR (Transparent)
# Back Surface: 50/50 Splitter Mirror
# Side Surface: Absorber (approximated as SimpleCoating(0.0, 0.0))

front_c = SimpleCoating(0.005)       # AR
back_c = SimpleCoating(0.5, 0.5)     # 50/50
side_c = SimpleCoating(0.0, 0.0)     # Absorber

# The Lens constructor automatically handles the spatial logic!
spatial_lens = Lens(
    SphericalSurface(50e-3, 25e-3),
    SphericalSurface(-50e-3, 25e-3),
    5e-3,
    x -> 1.5,
    front_c,
    back_c,
    side_c
)

println("Created lens with spatial coatings.")
println("Wrapper Type: $(typeof(spatial_lens.coating))")

# Let's inspect the coating at specific points
# We can use get_coating_at with local coordinates.
# Front Apex (Local 0,0,0 approx)
c_front = get_coating_at(
    spatial_lens.coating, spatial_lens, Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 1.0))
println("Coating at Front (Local 0,0,0): $(c_front)")

# Back Apex (Local 0, thickness, 0)
thick = 5e-3
c_back = get_coating_at(spatial_lens.coating, spatial_lens,
    Point3(0.0, thick + 1e-3, 0.0), Point3(0.0, 0.0, 1.0))
println("Coating at Back (Local y>thick): $(c_back)")

println("Coating at Back (Local y>thick): $(c_back)")

# ------------------------------------------------------------------------------
# 4. Meniscus Lens (Generic Spatial Logic)
# ------------------------------------------------------------------------------
println("\n--- 4. Meniscus Lens (Generic Logic) ---")
# Demonstrate that the new logic works for meniscus lenses where simple Y-split fails.
# Front: Convex (R=50mm), Back: Convex (R=100mm) -> Positive Meniscus
r1 = 50e-3
r2 = 100e-3
t = 5e-3
meniscus = Lens(
    SphericalSurface(r1, 25e-3),
    SphericalSurface(r2, 25e-3),
    t,
    x -> 1.5,
    SimpleCoating(0.1), # Front
    SimpleCoating(0.9)  # Back
)
# Check coating at Front Apex
c_meniscus_front = get_coating_at(
    meniscus.coating, meniscus, Point3(0.0, 0.0, 0.0), Point3(0.0, -1.0, 0.0))
println("Meniscus Front (R=$(r1*1e3)mm): $c_meniscus_front")

# Check coating at Back Apex
c_meniscus_back = get_coating_at(
    meniscus.coating, meniscus, Point3(0.0, t, 0.0), Point3(0.0, 1.0, 0.0))
println("Meniscus Back (R=$(r2*1e3)mm): $c_meniscus_back")

println("\nDone.")
