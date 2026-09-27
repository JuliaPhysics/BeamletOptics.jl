# Collimated ray bundle through a plano-convex singlet; find best focus by scanning a Detector.
# Run headless: julia --project=<env with BeamletOptics> 01_singlet_spot_diagram.jl
using BeamletOptics

const BMO = BeamletOptics
const mm = 1e-3
const nm = 1e-9

λ = 532nm
N_BK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945,
                          0.00600069867, 0.0200179144, 103.560653)

# Thorlabs LA1805-like plano-convex lens, curved side facing the source.
# Front vertex at y = 0, optical axis +y, ISO 10110 radius sign (center of curvature at +y -> R > 0).
r1, r2, ct, d = 15.5mm, Inf, 8.6mm, 25.4mm
lens = SphericalLens(r1, r2, ct, d, N_BK7)

n = N_BK7(λ)
f = BMO.lensmakers_eq(r1, r2, n)               # returns f (not 1/f), thin-lens estimate
bfl = f - ct / n                               # back focal length of a plano-convex lens
det = Detector(10mm)

system = System([lens, det])
source = CollimatedSource([0, -20mm, 0], [0, 1, 0], 6mm, λ; num_rings = 10, num_rays = 2000)

function rms_radius(spots)
    xs = first.(spots); zs = last.(spots)
    cx = sum(xs) / length(xs); cz = sum(zs) / length(zs)
    return sqrt(sum(@. (xs - cx)^2 + (zs - cz)^2) / length(spots))
end

best = (y = NaN, rms = Inf)
for y in range(ct + bfl - 2mm, ct + bfl + 1mm, length = 31)
    translate_to3d!(det, [0, y, 0])            # absolute move (translate3d! is relative)
    empty!(det)                                # detectors accumulate hits -> reset every run
    solve_system!(system, source)
    spots = spot_diagram(det)                  # Vector{Point2}: local (x, z) hit points in m
    r = rms_radius(spots)
    r < best.rms && (global best = (y = y, rms = r))
end

println("thin-lens f        = ", round(f / mm; digits = 3), " mm")
println("paraxial BFL       = ", round(bfl / mm; digits = 3), " mm")
println("best focus (vertex)= ", round((best.y - ct) / mm; digits = 3), " mm behind back vertex")
println("min RMS spot radius= ", round(best.rms / 1e-6; digits = 2), " µm")
