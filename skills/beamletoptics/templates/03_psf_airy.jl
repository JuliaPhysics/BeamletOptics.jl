# Diffraction-limited PSF of a singlet: coherent ray summation on a Detector.
using BeamletOptics

const mm = 1e-3
const µm = 1e-6

λ = 1000e-9
lens = SphericalLens(100mm, Inf, 1mm, 25.4mm, λ -> 1.5)

det = Detector(10mm)
translate3d!(det, [0, 200mm + 0.13mm, 0])     # near best focus

system = System([lens, det])
# UniformDiscSource: equal-area sampling, the right choice for PSF/intensity work
source = UniformDiscSource([0, -10mm, 0], [0, 1, 0], 15mm, λ; num_rays = 5000)
solve_system!(system, source)

x, z, I = intensity(det; n = 301, crop_factor = 10)   # grid around the spot, I in W/m²
imax = argmax(I)
println("peak at x = ", round(x[imax[1]] / µm; digits = 2), " µm, z = ",
        round(z[imax[2]] / µm; digits = 2), " µm")
println("Airy radius (1.22 λ N) ≈ ", round(1.22 * λ * 200mm / 15mm / µm; digits = 2), " µm")
