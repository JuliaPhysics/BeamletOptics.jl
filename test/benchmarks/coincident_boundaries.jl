# Runtime of tracing through components with touching elements: a cemented doublet, a plate and a cube
# beamsplitter, for Beam, GaussianBeamlet and AstigmaticGaussianBeamlet.
# Not part of the test suite. Run with
#
#   julia --project=test test/benchmarks/coincident_boundaries.jl
#
# Each case is solved once for warm-up, then `N` times, the minimum is reported (in µs per solve_system! call).
# Compare the numbers before and after a change of the tracing of touching elements.

using BeamletOptics

const mm = 1e-3
const λ = 1e-6
const N = 200

doublet() = begin
    dl = SphericalDoubletLens(50mm, -40mm, -100mm, 4mm, 2mm, 20mm, _ -> 1.5, _ -> 1.7)
    translate3d!(dl, [0, 30mm, 0])
    dl
end
plate() = begin
    pbs = RectangularPlateBeamsplitter(36mm, 25mm, 3mm, _ -> 1.5; reflectance = 0.3)
    zrotate3d!(pbs, deg2rad(30))
    pbs
end
cube() = CubeBeamsplitter(20mm, 1.5; reflectance = 0.3)

components = (
    "doublet" => (doublet, [2mm, -20mm, 0.0]),
    "plate beamsplitter" => (plate, [0.0, -20mm, 0.0]),
    "cube beamsplitter" => (cube, [3mm, -30mm, 0.0]))

sources = (
    "Beam" => (p) -> Beam(p, [0.0, 1, 0], λ),
    "GaussianBeamlet" => (p) -> GaussianBeamlet(p, [0.0, 1, 0], λ, 0.3mm),
    "AstigmaticGaussianBeamlet" => (p) -> AstigmaticGaussianBeamlet(p, [0.0, 1, 0], λ, 0.3mm))

function min_time(system, beam, n)
    solve_system!(system, beam) # warm-up, compiles
    t = Inf
    for _ in 1:n
        t = min(t, @elapsed solve_system!(system, beam))
    end
    return t
end

println("Julia ", VERSION, ", threads: ", Threads.nthreads(), ", ", Sys.cpu_info()[1].model, ", N = ", N)
println()
println(rpad("component", 22), rpad("source", 28), rpad("System [µs]", 14), "StaticSystem [µs]")
# 50 lenses around the component, off the beam path: the cost of the other objects of a system
function crowd()
    return [begin
                l = SphericalLens(30mm, -30mm, 5mm, 15mm, _ -> 1.5)
                φ = 2π * k / 50
                translate3d!(l, [80mm * cos(φ), 80mm * sin(φ), 20mm])
                l
            end
            for k in 1:50]
end

for crowded in (false, true), (cname, (make, p)) in components, (sname, source) in sources
    obj = make()
    objs = crowded ? [obj; crowd()] : [obj]
    cname = crowded ? cname * " + 50" : cname
    t_sys = min_time(System(objs), source(p), N)
    t_static = min_time(StaticSystem(objs), source(p), N)
    println(rpad(cname, 22), rpad(sname, 28), rpad(round(1e6 * t_sys; digits = 1), 14), round(1e6 * t_static; digits = 1))
end
