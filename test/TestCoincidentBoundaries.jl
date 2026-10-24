module TestCoincidentBoundaries

using BeamletOptics
using LinearAlgebra
using Test

const BMO = BeamletOptics

const mm = 1e-3

#=
Test matrix for touching elements: cemented lenses (DoubletLens, TripletLens), plate and cube
beamsplitters and bodies that merely touch in a System. At a shared boundary two or more surfaces
coincide, and the solver has to decide which interaction happens (today: per-component `Hint`s).

Every expectation is derived without BMO: the vector form of Snell's law and the law of reflection
in `follow`, an exact meridional trace through spherical surfaces (`trace_surfaces`), closed-form
results (plate offset, Fresnel power transmission, ABCD law of Gaussian beams) and physical
constructions (equal glass = no interface, reversibility of light paths, covariance under rigid
motion). Touching elements need no `Hint`: the solver collects the surfaces at a hit point
(`coincident_hits`) and derives the media on both sides from them (`transition`).
=#

#=
Helpers
=#

λ0 = 587.6e-9

"""Depth-first list of all beams (or beamlets) of a tree, the root first."""
all_beams(b) = vcat([b], reduce(vcat, all_beams.(BMO.children(b)); init = typeof(b)[]))
leaves(b) = filter(x -> isempty(BMO.children(x)), all_beams(b))

# the beam that carries the geometric path
chief_of(b::Beam) = b
chief_of(g::GaussianBeamlet) = g.chief
chief_of(a::AstigmaticGaussianBeamlet) = a.c
# all component beams of a beamlet
components(g::GaussianBeamlet) = (g.chief, g.waist, g.divergence)
components(a::AstigmaticGaussianBeamlet) = BMO._component_beams(a)

"""Position, direction and refractive index of every ray of every beam of the tree as one flat vector."""
function flatten(b)
    out = Float64[]
    for x in all_beams(b), r in rays(chief_of(x))
        append!(out, position(r), direction(r), BMO.refractive_index(r))
    end
    return out
end

"""Hit shape of every ray of a beam, `nothing` for a ray without intersection."""
shape_sequence(beam::Beam) = [isnothing(BMO.intersection(r)) ? nothing : BMO.shape(BMO.intersection(r))
                              for r in rays(beam)]

same_sequence(a, b) = length(a) == length(b) && all(x === y for (x, y) in zip(a, b))

# fresh beam from `pos` along `dir`
beam_at(pos, dir; λ = λ0) = Beam(collect(Float64, pos), collect(Float64, dir), λ)

function solve(objects, beam; kwargs...)
    solve_system!(System(objects), beam; kwargs...)
    return beam
end

function end_point(beam)
    r = last(rays(beam))
    return position(r), direction(r)
end

#=
Reference 1: vector Snell/reflection in `follow`. Surfaces are events along a ray; a surface is a plane or a sphere.
=#

"""
A surface to cross: `t(p, d)` is the ray parameter of the hit, `normal(point)` a normal of the surface at
the hit, `n` the refractive index behind it. `mirror` surfaces reflect, all others refract or, beyond the
critical angle, reflect (total internal reflection).
"""
struct Surface
    t::Function
    normal::Function
    n::Float64
    mirror::Bool
end

function plane(p0, N, n; mirror = false)
    N = normalize(collect(Float64, N))
    p0 = collect(Float64, p0)
    return Surface((p, d) -> dot(p0 - p, N) / dot(d, N), _ -> N, n, mirror)
end

"""Sphere surface with the vertex at height `yv` on the y axis and the radius of curvature `R` (BMO sign convention); `Inf` is a plane."""
function sphere(yv, R, n)
    isinf(R) && return plane([0, yv, 0], [0, 1, 0], n)
    c = [0.0, yv + R, 0.0]
    function t(p, d)
        oc = p - c
        b = dot(oc, d)
        disc = b^2 - (dot(oc, oc) - R^2)
        disc < 0 && return NaN
        ts = (-b - sqrt(disc), -b + sqrt(disc))
        ys = [abs((p + s * d)[2] - yv) for s in ts]
        return ts[argmin(ys)]
    end
    return Surface(t, pt -> normalize(pt - c), n, false)
end

"""
    follow(surfaces, p, d; n0 = 1.0)

Follows the ray from `p` along `d` through the `surfaces` in the given order. Returns the hits
`(pos, dir, n, tir)` behind each surface and the optical path length `opl` up to the last surface.
"""
function follow(surfaces, p, d; n0 = 1.0)
    p = collect(Float64, p)
    d = normalize(collect(Float64, d))
    n = n0
    hits = NamedTuple[]
    opl = 0.0
    for s in surfaces
        t = s.t(p, d)
        p = p + t * d
        opl += n * t
        N = s.normal(p)
        dot(N, d) > 0 && (N = -N)
        tir = false
        if s.mirror
            d = d - 2dot(d, N) * N
        else
            μ = n / s.n
            c = -dot(N, d)
            k = 1 - μ^2 * (1 - c^2)
            if k < 0
                tir = true
                d = d - 2dot(d, N) * N
            else
                d = normalize(μ * d + (μ * c - sqrt(k)) * N)
                n = s.n
            end
        end
        push!(hits, (; pos = copy(p), dir = copy(d), n, tir))
    end
    return (; hits, opl)
end

"""Surfaces of a stack of spherical elements, `vertices[i]` is the vertex of surface `i` on the y axis (forward)."""
lens_surfaces(vertices, radii, ns) = [sphere(y, R, n) for (y, R, n) in zip(vertices, radii, ns)]

"""The same surfaces seen from the other side, for a ray along -y (mirror image y -> -y)."""
function reversed_lens_surfaces(vertices, radii, ns; n0 = 1.0)
    after = [reverse(ns[1:(end - 1)]); n0]
    return [sphere(-y, -R, n) for (y, R, n) in zip(reverse(vertices), reverse(radii), after)]
end

# follow a ray along -y through `reversed_lens_surfaces`, result in the unmirrored frame
function follow_reversed(vertices, radii, ns, x0, y0)
    ref = follow(reversed_lens_surfaces(vertices, radii, ns), [x0, -y0, 0], [0, 1, 0])
    hits = [(; h..., pos = h.pos .* [1, -1, 1], dir = h.dir .* [1, -1, 1]) for h in ref.hits]
    return (; hits, opl = ref.opl)
end

"""
Largest deviation `(position [m], direction, refractive index)` between the rays of a single-leaf `beam` and a
`follow` result: ray `k + 1` of the beam starts at hit `k`. A different number of rays is `Inf`.
"""
function path_error(beam, ref)
    r = rays(beam)
    length(r) == length(ref.hits) + 1 || return (Inf, Inf, Inf)
    ep = ed = en = 0.0
    for (k, h) in enumerate(ref.hits)
        ray = r[k + 1]
        ep = max(ep, norm(position(ray) - h.pos))
        ed = max(ed, norm(direction(ray) - h.dir))
        en = max(en, abs(BMO.refractive_index(ray) - h.n))
    end
    return (ep, ed, en)
end

const TOL_POS = 1e-8
const TOL_DIR = 1e-7
agrees(err; dir = TOL_DIR) = err[1] < TOL_POS && err[2] < dir && err[3] < 1e-12

#=
Reference 2: the components under test
=#

# crown/flint doublet, positive ROC = center at +y
const DOUBLET = (r = (50mm, -40mm, -100mm), l = (4mm, 2mm), d = 20mm)
# crown/flint/crown triplet
const TRIPLET = (r = (40mm, -35mm, 60mm, -120mm), l = (3mm, 1.5mm, 2.5mm), d = 20mm)

function doublet(n1, n2; y0 = 0.0, g = DOUBLET)
    dl = SphericalDoubletLens(g.r..., g.l..., g.d, λ -> n1, λ -> n2)
    translate3d!(dl, [0, y0, 0])
    return dl
end
doublet_vertices(y0; g = DOUBLET) = [y0, y0 + g.l[1], y0 + sum(g.l)]
doublet_surfaces(y0, n1, n2; g = DOUBLET) = (doublet_vertices(y0; g), collect(g.r), [n1, n2, 1.0])

function triplet(n1, n2, n3; y0 = 0.0, g = TRIPLET)
    tl = SphericalTripletLens(g.r..., g.l..., g.d, λ -> n1, λ -> n2, λ -> n3)
    translate3d!(tl, [0, y0, 0])
    return tl
end
triplet_vertices(y0; g = TRIPLET) = y0 .+ [0.0, g.l[1], g.l[1] + g.l[2], sum(g.l)]
triplet_surfaces(y0, n1, n2, n3; g = TRIPLET) = (triplet_vertices(y0; g), collect(g.r), [n1, n2, n3, 1.0])

const HEIGHTS = (0.0, 1mm, -2mm, 3.5mm)

#=
Tests
=#

@testset "Coincident boundaries" begin
    @testset "1. Cemented SDF parts: doublet and triplet" begin
        y0 = 30mm
        ystart = -20mm
        @testset "Equal glass: a doublet is the equivalent singlet" begin
            n = 1.6
            dl = doublet(n, n; y0)
            single = SphericalLens(DOUBLET.r[1], DOUBLET.r[3], sum(DOUBLET.l), DOUBLET.d, λ -> n)
            translate3d!(single, [0, y0, 0])
            for h in HEIGHTS
                bd = solve([dl], beam_at([h, ystart, 0], [0, 1, 0]))
                bs = solve([single], beam_at([h, ystart, 0], [0, 1, 0]))
                pd, dd = end_point(bd)
                ps, ds = end_point(bs)
                @test pd ≈ ps atol = TOL_POS
                @test dd ≈ ds atol = TOL_DIR
                @test BMO.optical_path_length(bd) ≈ BMO.optical_path_length(bs) atol = TOL_POS
                # the cemented surface is no boundary: the doublet has one ray more, at the same index
                @test length(rays(bd)) == length(rays(bs)) + 1
                @test BMO.refractive_index.(rays(bd)) == [1.0, n, n, 1.0]
            end
        end

        @testset "Equal glass: a triplet is the equivalent singlet" begin
            n = 1.55
            tl = triplet(n, n, n; y0)
            single = SphericalLens(TRIPLET.r[1], TRIPLET.r[4], sum(TRIPLET.l), TRIPLET.d, λ -> n)
            translate3d!(single, [0, y0, 0])
            for h in HEIGHTS
                bt = solve([tl], beam_at([h, ystart, 0], [0, 1, 0]))
                bs = solve([single], beam_at([h, ystart, 0], [0, 1, 0]))
                @test end_point(bt)[1] ≈ end_point(bs)[1] atol = TOL_POS
                @test end_point(bt)[2] ≈ end_point(bs)[2] atol = TOL_DIR
                @test BMO.optical_path_length(bt) ≈ BMO.optical_path_length(bs) atol = TOL_POS
                @test BMO.refractive_index.(rays(bt)) == [1.0, n, n, n, 1.0]
            end
        end

        @testset "Triplet with two equal front glasses is the doublet" begin
            n1, n3 = 1.5, 1.72
            tl = triplet(n1, n1, n3; y0)
            g = (r = (TRIPLET.r[1], TRIPLET.r[3], TRIPLET.r[4]), l = (TRIPLET.l[1] + TRIPLET.l[2], TRIPLET.l[3]), d = TRIPLET.d)
            dl = doublet(n1, n3; y0, g)
            for h in HEIGHTS
                bt = solve([tl], beam_at([h, ystart, 0], [0, 1, 0]))
                bd = solve([dl], beam_at([h, ystart, 0], [0, 1, 0]))
                @test end_point(bt)[1] ≈ end_point(bd)[1] atol = TOL_POS
                @test end_point(bt)[2] ≈ end_point(bd)[2] atol = TOL_DIR
            end
        end

        @testset "Tilted plano doublet: Snell at the cemented surface by hand" begin
            n1, n2 = 1.5, 1.72
            l1, l2 = 4mm, 3mm
            dl = SphericalDoubletLens(Inf, Inf, Inf, l1, l2, 60mm, λ -> n1, λ -> n2)
            α = deg2rad(25)
            zrotate3d!(dl, α)
            beam = solve([dl], beam_at([0, -10mm, 0], [0, 1, 0]))
            # the lens axis u and the in-plane unit vector v, d = cos(α) u + sin(α) v
            u = [-sin(α), cos(α), 0]
            v = [cos(α), sin(α), 0]
            θ1 = asin(sin(α) / n1)
            θ2 = asin(sin(α) / n2)
            r = rays(beam)
            @test length(r) == 4
            # hit points: along u the layer thickness, along v the lateral walk
            p1 = zeros(3)
            p2 = l1 * u + l1 * tan(θ1) * v
            p3 = p2 + l2 * u + l2 * tan(θ2) * v
            @test position(r[2]) ≈ p1 atol = TOL_POS
            @test position(r[3]) ≈ p2 atol = TOL_POS
            @test position(r[4]) ≈ p3 atol = TOL_POS
            # directions behind the surfaces
            @test direction(r[2]) ≈ cos(θ1) * u + sin(θ1) * v atol = TOL_DIR
            @test direction(r[3]) ≈ cos(θ2) * u + sin(θ2) * v atol = TOL_DIR
            @test direction(r[4]) ≈ [0, 1, 0] atol = TOL_DIR
            @test BMO.refractive_index.(r) == [1.0, n1, n2, 1.0]
            # optical path length: 10 mm of air, then n l / cos(θ) per layer
            @test BMO.optical_path_length(beam) ≈ 10mm + n1 * l1 / cos(θ1) + n2 * l2 / cos(θ2) atol = TOL_POS
        end

        @testset "Curved doublet against the exact trace" begin
            n1, n2 = 1.5, 1.72
            dl = doublet(n1, n2; y0)
            for h in HEIGHTS
                beam = solve([dl], beam_at([h, ystart, 0], [0, 1, 0]))
                ref = follow(lens_surfaces(doublet_surfaces(y0, n1, n2)...), [h, ystart, 0], [0, 1, 0])
                @test agrees(path_error(beam, ref))
                @test BMO.optical_path_length(beam) ≈ ref.opl atol = TOL_POS
            end
        end

        @testset "Triplet against the exact trace, three different glasses" begin
            n = (1.5, 1.74, 1.62)
            tl = triplet(n...; y0)
            for h in HEIGHTS
                beam = solve([tl], beam_at([h, ystart, 0], [0, 1, 0]))
                ref = follow(lens_surfaces(triplet_surfaces(y0, n...)...), [h, ystart, 0], [0, 1, 0])
                @test agrees(path_error(beam, ref))
                @test BMO.optical_path_length(beam) ≈ ref.opl atol = TOL_POS
            end
        end

        @testset "Skew ray through the doublet" begin
            n1, n2 = 1.5, 1.72
            dl = doublet(n1, n2; y0)
            d = normalize([0.05, 1, -0.03])
            p = [-3mm, ystart, 1.5mm]
            beam = solve([dl], beam_at(p, d))
            ref = follow(lens_surfaces(doublet_surfaces(y0, n1, n2)...), p, d)
            @test agrees(path_error(beam, ref))
        end
    end

    @testset "2. Cemented mesh parts" begin
        # two mesh-based refractive plates (CuboidMesh prisms), the first face of each at its y coordinate
        nA, nB = 1.7, 1.4
        tA, tB = 5mm, 3mm
        θ = deg2rad(30)
        d = [sin(θ), cos(θ), 0]
        p = [0, -10mm, 0]
        function plates(gap)
            A = RectangularCompensatorPlate(40mm, 40mm, tA, λ -> nA)
            B = RectangularCompensatorPlate(40mm, 40mm, tB, λ -> nB)
            translate3d!(B, [0, tA + gap, 0])
            return A, B
        end
        function plate_surfaces(gap)
            return [plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], 1.0 * (gap == 0 ? nB : 1.0)),
                plane([0, tA + gap, 0], [0, 1, 0], nB), plane([0, tA + gap + tB, 0], [0, 1, 0], 1.0)]
        end

        @testset "Single plate" begin
            A, _ = plates(0.0)
            beam = solve([A], beam_at(p, d))
            ref = follow([plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], 1.0)], p, d)
            @test agrees(path_error(beam, ref))
        end

        @testset "Separated by an air gap: refraction into air, then into the second plate" begin
            A, B = plates(1mm)
            beam = solve([A, B], beam_at(p, d))
            ref = follow(plate_surfaces(1mm), p, d)
            @test agrees(path_error(beam, ref))
            @test BMO.refractive_index.(rays(beam)) == [1.0, nA, 1.0, nB, 1.0]
            # the exit is parallel to the entry (plane parallel layers)
            @test direction(last(rays(beam))) ≈ d atol = TOL_DIR
        end

        @testset "Touching: Snell from the first into the second glass" begin
            A, B = plates(0.0)
            beam = solve([A, B], beam_at(p, d))
            # first glass into second glass, no air in between
            ref = follow([plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], nB),
                    plane([0, tA + tB, 0], [0, 1, 0], 1.0)], p, d)
            err = path_error(beam, ref)
            # the plates are separate objects: the rear face of A and the front face of B are one boundary
            @test agrees(err)
            @test BMO.refractive_index.(rays(beam)) == [1.0, nA, nB, 1.0]
            @test direction(rays(beam)[3]) ≈ ref.hits[2].dir atol = TOL_DIR
        end

        @testset "Three touching plates" begin
            nC, tC = 1.55, 2mm
            A = RectangularCompensatorPlate(40mm, 40mm, tA, λ -> nA)
            B = RectangularCompensatorPlate(40mm, 40mm, tB, λ -> nB)
            C = RectangularCompensatorPlate(40mm, 40mm, tC, λ -> nC)
            translate3d!(B, [0, tA, 0])
            translate3d!(C, [0, tA + tB, 0])
            beam = solve([A, B, C], beam_at(p, d))
            ref = follow([plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], nB),
                    plane([0, tA + tB, 0], [0, 1, 0], nC), plane([0, tA + tB + tC, 0], [0, 1, 0], 1.0)], p, d)
            # every touching face pair is one boundary
            @test agrees(path_error(beam, ref))
        end
    end

    @testset "3. Beamsplitters: paths, split and leaves" begin
        R = 0.3
        n = 1.5
        @testset "Plate beamsplitter, coating side first" begin
            t = 3mm
            α = deg2rad(40)
            pbs = RectangularPlateBeamsplitter(36mm, 25mm, t, λ -> n; reflectance = R)
            zrotate3d!(pbs, α)
            E0 = [1.0, 0, 0]
            beam = Beam([0.0, -20mm, 0], [0.0, 1, 0], λ0, E0)
            solve_system!(System([pbs]), beam)
            # coating normal (faces the beam) and the substrate behind it
            N = [sin(α), -cos(α), 0]
            @test length(all_beams(beam)) == 3
            @test length(leaves(beam)) == 2
            tr, rf = BMO.children(beam)
            # reflected leaf: one ray, law of reflection
            d0 = [0.0, 1, 0]
            @test length(rays(rf)) == 1
            @test position(first(rays(rf))) ≈ zeros(3) atol = TOL_POS
            @test direction(first(rays(rf))) ≈ d0 - 2dot(d0, N) * N atol = TOL_DIR
            # transmitted leaf: refracted into the substrate, then out of its back face
            ref = follow([plane([0, 0, 0], N, n), plane(-t * N, N, 1.0)], [0, -20mm, 0], d0)
            @test length(rays(tr)) == 2
            @test position(rays(tr)[1]) ≈ zeros(3) atol = TOL_POS
            @test direction(rays(tr)[1]) ≈ ref.hits[1].dir atol = TOL_DIR
            @test position(rays(tr)[2]) ≈ ref.hits[2].pos atol = TOL_POS
            @test direction(rays(tr)[2]) ≈ d0 atol = TOL_DIR
            @test BMO.refractive_index.(rays(tr)) == [n, 1.0]
            # textbook lateral displacement of a tilted plate: t sin(θ - θ') / cos(θ')
            θ = α
            θ′ = asin(sin(θ) / n)
            @test abs(position(rays(tr)[2])[1]) ≈ t * sin(θ - θ′) / cos(θ′) atol = TOL_POS
            # power: the reflection takes R, the transmission 1 - R. Power ~ n cos(θ) |E|²
            P0 = norm(E0)^2
            @test sum(abs2, BMO.polarization(first(rays(rf)))) ≈ R * P0 atol = 1e-12
            Pt = n * cos(θ′) * sum(abs2, BMO.polarization(rays(tr)[1]))
            @test Pt ≈ (1 - R) * cos(θ) * P0 atol = 1e-12
        end

        @testset "Cube beamsplitter, tilted" begin
            L = 20mm
            ρ = deg2rad(20)
            cbs = CubeBeamsplitter(L, n; reflectance = R)
            zrotate3d!(cbs, ρ)
            Rz = [cos(ρ) -sin(ρ) 0; sin(ρ) cos(ρ) 0; 0 0 1]
            p0 = [2mm, -30mm, 0.0]
            d0 = [0.0, 1, 0]
            beam = solve([cbs], beam_at(p0, d0))
            # faces of the unrotated cube centered at the origin and the hypotenuse x + y = 0
            Nh = Rz * normalize([1, 1, 0])
            Nbot, Ntop = Rz * [0, -1, 0], Rz * [0, 1, 0]
            Nleft = Rz * [-1, 0, 0]
            @test length(leaves(beam)) == 2
            tr, rf = BMO.children(beam)
            entry = plane(L / 2 * Nbot, Nbot, n)
            # transmitted: entry face, coating (same glass, no change), exit face
            ref_t = follow([entry, plane(zeros(3), Nh, n), plane(L / 2 * Ntop, Ntop, 1.0)], p0, d0)
            # reflected: entry face, mirror at the coating, left face
            ref_r = follow([entry, plane(zeros(3), Nh, n; mirror = true), plane(L / 2 * Nleft, Nleft, 1.0)], p0, d0)
            @test agrees(path_error(beam, (; hits = ref_t.hits[1:1])))
            @test length(rays(beam)) == 2
            # the first ray of each leaf starts at the coating, the last one at the exit face
            @test position(first(rays(tr))) ≈ ref_t.hits[2].pos atol = TOL_POS
            @test position(last(rays(tr))) ≈ ref_t.hits[3].pos atol = TOL_POS
            @test direction(last(rays(tr))) ≈ d0 atol = TOL_DIR
            @test position(first(rays(rf))) ≈ ref_r.hits[2].pos atol = TOL_POS
            @test direction(first(rays(rf))) ≈ ref_r.hits[2].dir atol = TOL_DIR
            @test position(last(rays(rf))) ≈ ref_r.hits[3].pos atol = TOL_POS
            @test direction(last(rays(rf))) ≈ ref_r.hits[3].dir atol = TOL_DIR
            @test length(rays(tr)) == 2
            @test length(rays(rf)) == 2
            @test BMO.refractive_index.(rays(beam)) == [1.0, n]
            @test BMO.refractive_index.(rays(tr)) == [n, 1.0]
            @test BMO.refractive_index.(rays(rf)) == [n, 1.0]
        end

        @testset "Cube beamsplitter, normal incidence: power split" begin
            L = 20mm
            cbs = CubeBeamsplitter(L, n; reflectance = R)
            E0 = [1.0, 0, 0]
            for x in (0.0, 3mm, -9.9mm, 9.9mm)
                beam = Beam([x, -30mm, 0], [0.0, 1, 0], λ0, E0)
                solve_system!(System([cbs]), beam)
                @test length(leaves(beam)) == 2
                tr, rf = BMO.children(beam)
                # hypotenuse x + y = 0: coating hit at (x, -x)
                @test position(first(rays(tr))) ≈ [x, -x, 0] atol = TOL_POS
                @test position(last(rays(tr))) ≈ [x, L / 2, 0] atol = TOL_POS
                @test position(last(rays(rf))) ≈ [-L / 2, -x, 0] atol = TOL_POS
                @test direction(last(rays(tr))) ≈ [0, 1, 0] atol = TOL_DIR
                @test direction(last(rays(rf))) ≈ [-1, 0, 0] atol = TOL_DIR
                # Fresnel at normal incidence: T = 4n/(1+n)² per uncoated face, coating splits R : 1 - R
                Tf = 4n / (1 + n)^2
                flux(ray) = BMO.refractive_index(ray) * sum(abs2, BMO.polarization(ray))
                @test flux(last(rays(tr))) ≈ Tf^2 * (1 - R) atol = 1e-12
                @test flux(last(rays(rf))) ≈ Tf^2 * R atol = 1e-12
            end
        end
    end

    @testset "4. Reversed passes" begin
        y0 = 30mm
        n1, n2 = 1.5, 1.72
        @testset "Doublet from the back against the exact trace" begin
            dl = doublet(n1, n2; y0)
            for h in HEIGHTS
                beam = solve([dl], beam_at([h, 80mm, 0], [0, -1, 0]))
                ref = follow_reversed(doublet_surfaces(y0, n1, n2)..., h, 80mm)
                @test agrees(path_error(beam, ref))
                @test BMO.optical_path_length(beam) ≈ ref.opl atol = TOL_POS
                @test BMO.refractive_index.(rays(beam)) == [1.0, 1.72, 1.5, 1.0]
            end
        end

        @testset "Triplet from the back against the exact trace" begin
            n = (1.5, 1.74, 1.62)
            tl = triplet(n...; y0)
            for h in HEIGHTS
                beam = solve([tl], beam_at([h, 80mm, 0], [0, -1, 0]))
                ref = follow_reversed(triplet_surfaces(y0, n...)..., h, 80mm)
                @test agrees(path_error(beam, ref))
            end
        end

        @testset "Reversibility of the light path through doublet and triplet" begin
            for (name, obj) in (("doublet", doublet(n1, n2; y0)), ("triplet", triplet(1.5, 1.74, 1.62; y0)))
                fwd = solve([obj], beam_at([2mm, -20mm, 0], [0, 1, 0]))
                p, d = end_point(fwd)
                # send the light back from a point behind the exit
                back = solve([obj], beam_at(p + 40mm * d, -d))
                pb, db = end_point(back)
                @test db ≈ [0, -1, 0] atol = TOL_DIR
                # same lateral height at the entry plane, the path is the same ray reversed
                @test pb[1] ≈ 2mm atol = TOL_POS
                @test BMO.refractive_index.(rays(back)) == reverse(BMO.refractive_index.(rays(fwd)))
                air_in = norm(position(rays(fwd)[2]) - position(rays(fwd)[1]))
                @test BMO.optical_path_length(back) - 40mm ≈ BMO.optical_path_length(fwd) - air_in atol = TOL_POS
            end
        end

        @testset "Plate beamsplitter from the substrate side" begin
            t, n, R = 3mm, 1.5, 0.3
            α = deg2rad(40)
            pbs = RectangularPlateBeamsplitter(36mm, 25mm, t, λ -> n; reflectance = R)
            zrotate3d!(pbs, α)
            N = [sin(α), -cos(α), 0]
            p0 = [0.0, 20mm, 0.0]
            d0 = [0.0, -1, 0]
            beam = solve([pbs], beam_at(p0, d0))
            back_face = plane(-t * N, N, n)
            coating = plane(zeros(3), N, 1.0)
            # transmitted: back face in, coating out into air. Reflected: coating mirrors, back face out
            ref_in = follow([back_face], p0, d0)
            ref_t = follow([back_face, coating], p0, d0)
            ref_r = follow([back_face, plane(zeros(3), N, n; mirror = true), plane(-t * N, N, 1.0)], p0, d0)
            @test length(rays(beam)) == 2
            @test position(rays(beam)[2]) ≈ ref_in.hits[1].pos atol = TOL_POS
            @test direction(rays(beam)[2]) ≈ ref_in.hits[1].dir atol = TOL_DIR
            @test length(leaves(beam)) == 2
            tr, rf = BMO.children(beam)
            # transmitted: leaves into air, plane parallel, so along the incident direction
            @test length(rays(tr)) == 1
            @test position(first(rays(tr))) ≈ ref_t.hits[2].pos atol = TOL_POS
            @test direction(first(rays(tr))) ≈ ref_t.hits[2].dir atol = TOL_DIR
            @test BMO.refractive_index(first(rays(tr))) == 1.0
            # reflected: inside, then out of the back face
            @test length(rays(rf)) == 2
            @test position(rays(rf)[1]) ≈ ref_r.hits[2].pos atol = TOL_POS
            @test direction(rays(rf)[1]) ≈ ref_r.hits[2].dir atol = TOL_DIR
            @test position(rays(rf)[2]) ≈ ref_r.hits[3].pos atol = TOL_POS
            @test direction(rays(rf)[2]) ≈ ref_r.hits[3].dir atol = TOL_DIR
        end

        @testset "Cube beamsplitter from the other faces" begin
            L, n, R = 20mm, 1.5, 0.3
            cbs = CubeBeamsplitter(L, n; reflectance = R)
            x = 3mm
            # from the top along -y: reflection to +x, transmission to -y
            beam = solve([cbs], beam_at([x, 30mm, 0], [0, -1, 0]))
            tr, rf = BMO.children(beam)
            @test position(last(rays(tr))) ≈ [x, -L / 2, 0] atol = TOL_POS
            @test direction(last(rays(tr))) ≈ [0, -1, 0] atol = TOL_DIR
            @test position(last(rays(rf))) ≈ [L / 2, -x, 0] atol = TOL_POS
            @test direction(last(rays(rf))) ≈ [1, 0, 0] atol = TOL_DIR
            # from the +x side along -x: reflection to +y, transmission to -x
            beam = solve([cbs], beam_at([30mm, x, 0], [-1, 0, 0]))
            tr, rf = BMO.children(beam)
            @test position(last(rays(tr))) ≈ [-L / 2, x, 0] atol = TOL_POS
            @test position(last(rays(rf))) ≈ [-x, L / 2, 0] atol = TOL_POS
            @test direction(last(rays(rf))) ≈ [0, 1, 0] atol = TOL_DIR
        end

        @testset "Reversibility of the cube beamsplitter" begin
            L, n = 20mm, 1.5
            cbs = CubeBeamsplitter(L, n; reflectance = 0.4)
            zrotate3d!(cbs, deg2rad(15))
            p0 = [1mm, -30mm, 0.0]
            fwd = solve([cbs], beam_at(p0, [0, 1, 0]))
            for leaf in leaves(fwd)
                p, d = end_point(leaf)
                back = solve([cbs], beam_at(p + 30mm * d, -d))
                # one of the leaves of the reversed beam returns to the source
                found = any(leaves(back)) do b
                    pb, db = end_point(b)
                    isapprox(db, [0, -1, 0]; atol = TOL_DIR) && isapprox(pb[1], p0[1]; atol = TOL_POS)
                end
                @test found
            end
        end
    end

    @testset "5. Grazing incidence and TIR at a cemented surface" begin
        # high-index crown against low-index second glass: sin(φ_c) = n2 / n1 = 0.7368, i.e. the rim ray at
        # h = R sin(φ_c) = 8.84 mm of a collimated beam is critical at the convex cemented surface
        n1, n2 = 1.9, 1.4
        R2 = -12mm
        g = (r = (Inf, R2, Inf), l = (6mm, 8mm), d = 20mm)
        dl = doublet(n1, n2; y0 = 0.0, g)
        surf = doublet_surfaces(0.0, n1, n2; g)
        φc_height = abs(R2) * n2 / n1
        @test 8.8mm < φc_height < 8.9mm

        @testset "Grazing transmission just below the critical angle" begin
            for h in (7mm, 8.5mm, 8.8mm)
                beam = solve([dl], beam_at([h, -10mm, 0], [0, 1, 0]))
                ref = follow(lens_surfaces(surf...), [h, -10mm, 0], [0, 1, 0])
                @test !any(x -> x.tir, ref.hits) # the reference transmits
                # the direction error is amplified by 1 / cos(θ') at an exit angle θ' close to 90°
                @test agrees(path_error(beam, ref); dir = 1e-6)
            end
        end

        @testset "Total internal reflection at the cemented surface" begin
            h = 9mm
            beam = solve([dl], beam_at([h, -10mm, 0], [0, 1, 0]))
            ref = follow(lens_surfaces(surf...), [h, -10mm, 0], [0, 1, 0])
            @test ref.hits[2].tir # the reference reflects at the cemented surface
            r = rays(beam)
            # the reflection: position and direction by the law of reflection
            @test position(r[3]) ≈ ref.hits[2].pos atol = TOL_POS
            @test direction(r[3]) ≈ ref.hits[2].dir atol = TOL_DIR
            # the reflected ray stays in the front lens (n1)
            @test BMO.refractive_index(r[3]) == n1
            # the reflected ray hits the cemented surface a second time: the second event of the reference
            ref2 = follow([sphere(6mm, R2, n2), plane([0, 0, 0], [0, 1, 0], 1.0)], ref.hits[2].pos, ref.hits[2].dir; n0 = n1)
            @test ref2.hits[1].tir
            # the position and direction of this second reflection and the exit through the first surface are right,
            # in spite of the wrong index of the segment between
            @test position(r[4]) ≈ ref2.hits[1].pos atol = TOL_POS
            @test direction(r[4]) ≈ ref2.hits[1].dir atol = TOL_DIR
            @test length(r) == 5
            @test position(r[5]) ≈ ref2.hits[2].pos atol = TOL_POS
            @test direction(r[5]) ≈ ref2.hits[2].dir atol = TOL_DIR
        end

        @testset "Reversed grazing pass: low to high index" begin
            for h in (7mm, 8.5mm, 9mm)
                beam = solve([dl], beam_at([h, 40mm, 0], [0, -1, 0]))
                ref = follow_reversed(surf..., h, 40mm)
                @test !any(x -> x.tir, ref.hits)
                @test agrees(path_error(beam, ref))
            end
        end

        @testset "Touching prisms of different glass (cube without coating)" begin
            L = 20mm
            cbs = CubeBeamsplitter(L, 1.5)
            Nh = normalize([1, 1, 0])
            x = 3mm
            p0 = [x, -30mm, 0.0]
            d0 = [0.0, 1, 0]
            function prism_pair(na, nb)
                front = Prism(BMO.shape(cbs.front), λ -> na)
                back = Prism(BMO.shape(cbs.back), λ -> nb)
                return [front, back]
            end
            for (na, nb, name) in ((2.0, 1.35, "TIR at the hypotenuse"), (1.6, 1.4, "refraction into lower index"),
                (1.4, 1.6, "refraction into higher index"))
                beam = solve(prism_pair(na, nb), beam_at(p0, d0))
                entry = plane([0, -L / 2, 0], [0, -1, 0], na)
                hyp = plane(zeros(3), Nh, nb)
                ref = follow([entry, hyp], p0, d0)
                tir = ref.hits[2].tir
                tir || push!(ref.hits, follow([entry, hyp, plane([0, L / 2, 0], [0, 1, 0], 1.0)], p0, d0).hits[3])
                tir && push!(ref.hits, follow([entry, plane(zeros(3), Nh, na; mirror = true),
                    plane([-L / 2, 0, 0], [-1, 0, 0], 1.0)], p0, d0).hits[3])
                err = path_error(beam, ref)
                if name == "TIR at the hypotenuse"
                    @test tir
                    # right by accident: the front prism leaves into air at 45° and is totally reflected as well
                    @test agrees(err)
                else
                    @test !tir
                    # without a Hint the front prism refracts into the touching back prism
                    @test agrees(err)
                end
            end
        end
    end

    @testset "6. Kinematics" begin
        y0 = 30mm
        n1, n2 = 1.5, 1.72
        h = 2mm
        @testset "Doublet parts moved apart and back" begin
            dl = doublet(n1, n2; y0)
            before = flatten(solve([dl], beam_at([h, -20mm, 0], [0, 1, 0])))
            translate3d!(dl.back, [0, 3mm, 0])
            translate3d!(dl.back, [0, -3mm, 0])
            @test flatten(solve([dl], beam_at([h, -20mm, 0], [0, 1, 0]))) ≈ before atol = TOL_POS
            rotate3d!(dl.back, [0.0, 0, 1], deg2rad(2))
            rotate3d!(dl.back, [0.0, 0, 1], -deg2rad(2))
            @test flatten(solve([dl], beam_at([h, -20mm, 0], [0, 1, 0]))) ≈ before atol = TOL_POS
            # sideways and back
            translate3d!(dl.front, [1mm, 0, 2mm])
            translate3d!(dl.front, [-1mm, 0, -2mm])
            @test flatten(solve([dl], beam_at([h, -20mm, 0], [0, 1, 0]))) ≈ before atol = TOL_POS
        end

        @testset "Doublet parts apart: the air gap is seen" begin
            gap = 3mm
            dl = doublet(n1, n2; y0)
            translate3d!(dl.back, [0, gap, 0])
            beam = solve([dl], beam_at([h, -20mm, 0], [0, 1, 0]))
            v, R, ns = doublet_surfaces(y0, n1, n2)
            # front surface in air, front rear surface into the gap (air), back front surface into n2, back surface
            surfaces = [sphere(v[1], R[1], n1), sphere(v[2], R[2], 1.0), sphere(v[2] + gap, R[2], n2),
                sphere(v[3] + gap, R[3], 1.0)]
            ref = follow(surfaces, [h, -20mm, 0], [0, 1, 0])
            # the Hint of the front lens points to the back lens, but the rear surface of the front lens comes first
            @test agrees(path_error(beam, ref))
            @test BMO.refractive_index.(rays(beam)) == [1.0, n1, 1.0, n2, 1.0]
        end

        @testset "Rigid motion of the whole component" begin
            # rotated and translated component with the beam moved along: all rays move rigidly
            Rm = BMO.rotate3d([0.3, 0.2, 1.0], deg2rad(37))
            shift = [10mm, -5mm, 7mm]
            n3 = 1.62
            builders = (
                "doublet" => () -> doublet(n1, n2; y0 = 0.0),
                "triplet" => () -> triplet(1.5, 1.74, n3; y0 = 0.0),
                "plate" => () -> RectangularPlateBeamsplitter(36mm, 25mm, 3mm, λ -> 1.5; reflectance = 0.3),
                "cube" => () -> CubeBeamsplitter(20mm, 1.5; reflectance = 0.3))
            starts = ("doublet" => [h, -20mm, 0.0], "triplet" => [h, -20mm, 0.0],
                "plate" => [h, -20mm, 0.0], "cube" => [h, -30mm, 0.0])
            for ((name, build), (_, p0)) in zip(builders, starts)
                obj = build()
                base = flatten(solve([obj], beam_at(p0, [0, 1, 0])))
                obj2 = build()
                rotate3d!(obj2, Rm, zeros(3))
                translate3d!(obj2, shift)
                moved = flatten(solve([obj2], beam_at(Rm * p0 + shift, Rm * [0, 1, 0])))
                @test length(moved) == length(base)
                expected = similar(base)
                for k in 1:7:length(base)
                    expected[k:(k + 2)] = Rm * base[k:(k + 2)] + shift
                    expected[(k + 3):(k + 5)] = Rm * base[(k + 3):(k + 5)]
                    expected[k + 6] = base[k + 6]
                end
                @test moved ≈ expected atol = 1e-7
            end
        end

        @testset "Plate beamsplitter: coating moved off the substrate" begin
            t, n, R = 3mm, 1.5, 0.3
            gap = 2mm
            pbs = RectangularPlateBeamsplitter(36mm, 25mm, t, λ -> n; reflectance = R)
            before = flatten(solve([pbs], beam_at([0, -20mm, 0], [0, 1, 0])))
            translate3d!(pbs.coating, [0, -gap, 0])
            beam = solve([pbs], beam_at([0, -20mm, 0], [0, 1, 0]))
            tr, rf = BMO.children(beam)
            # the coating is in air, 2 mm in front of the substrate: the transmitted light travels in air, enters, leaves
            @test position(first(rays(rf))) ≈ [0, -gap, 0] atol = TOL_POS
            @test direction(first(rays(rf))) ≈ [0, -1, 0] atol = TOL_DIR
            # the transmitted ray stays in air between the coating and the substrate
            @test BMO.refractive_index(first(rays(tr))) == 1.0
            @test BMO.refractive_index.(rays(tr)) == [1.0, n, 1.0]
            # and back
            translate3d!(pbs.coating, [0, gap, 0])
            @test flatten(solve([pbs], beam_at([0, -20mm, 0], [0, 1, 0]))) ≈ before atol = TOL_POS
        end

        @testset "Cube beamsplitter: front prism moved off the coating" begin
            L, n = 20mm, 1.5
            cbs = CubeBeamsplitter(L, n; reflectance = 0.3)
            p0 = [3mm, -30mm, 0.0]
            before = flatten(solve([cbs], beam_at(p0, [0, 1, 0])))
            v = normalize([1.0, 1, 0])
            gap = 1mm
            translate3d!(cbs.front, -gap * v)
            beam = solve([cbs], beam_at(p0, [0, 1, 0]))
            # with an air gap the light at 45° inside the n = 1.5 prism is totally reflected at its hypotenuse
            # (n sin 45° = 1.06 > 1): reflected towards -x, leaves through the left leg, no splitting.
            ref = follow([plane([0, -L / 2 - gap * v[2], 0], [0, -1, 0], n),
                    plane(-gap * v, v, n; mirror = true), plane([-L / 2 - gap * v[1], 0, 0], [-1, 0, 0], 1.0)],
                p0, [0, 1, 0])
            # the Hint of the front prism points to the coating, but the hypotenuse of the prism comes first
            @test isempty(BMO.children(beam))
            @test position(last(rays(beam))) ≈ ref.hits[3].pos atol = TOL_POS
            @test direction(last(rays(beam))) ≈ ref.hits[3].dir atol = TOL_DIR
            # and back
            translate3d!(cbs.front, gap * v)
            @test flatten(solve([cbs], beam_at(p0, [0, 1, 0]))) ≈ before atol = TOL_POS
        end
    end

    @testset "7. More than two coincident hits" begin
        @testset "Cube hypotenuse: front, back and coating coincide" begin
            L, n, R = 20mm, 1.5, 0.3
            cbs = CubeBeamsplitter(L, n; reflectance = R)
            # from the center to the corners where the hypotenuse meets the faces
            for x in (0.0, 5mm, -5mm, 9.99mm, -9.99mm), dir in ([0.0, 1, 0], [0.0, -1, 0])
                p0 = [x, -dir[2] * 30mm, 0.0]
                beam = solve([cbs], beam_at(p0, dir))
                @test length(leaves(beam)) == 2
                @test all(b -> length(rays(b)) ≤ 2, all_beams(beam))
                tr, rf = BMO.children(beam)
                y_hit = -x
                @test position(first(rays(tr))) ≈ [x, y_hit, 0] atol = TOL_POS
                @test position(first(rays(rf))) ≈ [x, y_hit, 0] atol = TOL_POS
                @test direction(last(rays(tr))) ≈ dir atol = TOL_DIR
                @test direction(last(rays(rf))) ≈ (dir[2] > 0 ? [-1, 0, 0] : [1, 0, 0]) atol = TOL_DIR
            end
        end

        @testset "Thin beamsplitter between two touching plates" begin
            nA, nB = 1.7, 1.4
            tA, tB = 5mm, 3mm
            A = RectangularCompensatorPlate(40mm, 40mm, tA, λ -> nA)
            B = RectangularCompensatorPlate(40mm, 40mm, tB, λ -> nB)
            translate3d!(B, [0, tA, 0])
            bs = ThinBeamsplitter(40mm, 40mm; reflectance = 0.5)
            translate3d!(bs, [0, tA, 0])
            θ = 0.5
            d = [sin(θ), cos(θ), 0]
            p = [0.0, -10mm, 0.0]
            beam = solve([A, bs, B], beam_at(p, d))
            # three coincident hits: the splitter between the plates splits the beam
            @test length(leaves(beam)) == 2
            # reference of the transmitted branch: A into B (Snell) and out of B, and of the reflected one: back out of A
            ref_t = follow([plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], nB),
                    plane([0, tA + tB, 0], [0, 1, 0], 1.0)], p, d)
            @test any(b -> isapprox(end_point(b)[2], ref_t.hits[3].dir; atol = TOL_DIR), leaves(beam))
            ref_r = follow([plane([0, 0, 0], [0, 1, 0], nA), plane([0, tA, 0], [0, 1, 0], nA; mirror = true),
                    plane([0, 0, 0], [0, 1, 0], 1.0)], p, d)
            @test any(b -> isapprox(end_point(b)[2], ref_r.hits[3].dir; atol = TOL_DIR), leaves(beam))
        end

        @testset "Triplet: two cemented surfaces in sequence" begin
            # both cemented surfaces are crossed in both directions with the same result as the exact trace
            n = (1.5, 1.74, 1.62)
            tl = triplet(n...; y0 = 0.0)
            for dir in (+1, -1), h in (0.0, 2mm, 3.5mm)
                v, R, ns = triplet_surfaces(0.0, n...)
                if dir == 1
                    beam = solve([tl], beam_at([h, -20mm, 0], [0, 1, 0]))
                    ref = follow(lens_surfaces(v, R, ns), [h, -20mm, 0], [0, 1, 0])
                else
                    beam = solve([tl], beam_at([h, 40mm, 0], [0, -1, 0]))
                    ref = follow_reversed(v, R, ns, h, 40mm)
                end
                @test agrees(path_error(beam, ref))
            end
        end
    end

    @testset "8. StaticSystem equals System" begin
        y0 = 30mm
        cases = (
            "doublet" => (doublet(1.5, 1.72; y0), [2mm, -20mm, 0.0]),
            "triplet" => (triplet(1.5, 1.74, 1.62; y0), [2mm, -20mm, 0.0]),
            "plate" => (RectangularPlateBeamsplitter(36mm, 25mm, 3mm, λ -> 1.5; reflectance = 0.3), [0.0, -20mm, 0.0]),
            "cube" => (CubeBeamsplitter(20mm, 1.5; reflectance = 0.3), [3mm, -30mm, 0.0]))
        for (name, (obj, p0)) in cases
            name == "plate" && zrotate3d!(obj, deg2rad(30))
            @testset "$name" begin
                b1 = beam_at(p0, [0, 1, 0])
                b2 = beam_at(p0, [0, 1, 0])
                solve_system!(System([obj]), b1)
                solve_system!(StaticSystem([obj]), b2)
                @test flatten(b1) == flatten(b2)
                @test length(all_beams(b1)) == length(all_beams(b2))
                # polarized ray
                e1 = Beam(p0, [0.0, 1, 0], λ0, [1.0, 0, 0])
                e2 = Beam(p0, [0.0, 1, 0], λ0, [1.0, 0, 0])
                solve_system!(System([obj]), e1)
                solve_system!(StaticSystem([obj]), e2)
                pol(b) = reduce(vcat, [collect(BMO.polarization(r)) for x in all_beams(b) for r in rays(x)])
                @test pol(e1) == pol(e2)
            end
        end
        @testset "Gaussian beamlet through the doublet" begin
            dl = doublet(1.5, 1.72; y0)
            g1 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.3mm)
            g2 = GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], 1e-6, 0.3mm)
            solve_system!(System([dl]), g1)
            solve_system!(StaticSystem([dl]), g2)
            @test flatten(g1) == flatten(g2)
            @test all(first(BMO.gauss_parameters(g1, z)) == first(BMO.gauss_parameters(g2, z)) for z in (50mm, 80mm))
        end
    end

    @testset "9. Beamlets and polarized rays" begin
        λ = 1e-6
        w0 = 0.3mm
        q0 = im * π * w0^2 / λ
        P(d) = [1 d; 0 1]
        S(R, n1, n2) = [1 0; (n1 - n2) / (n2 * R) n1 / n2]
        prop(M, q) = (M[1, 1] * q + M[1, 2]) / (M[2, 1] * q + M[2, 2])
        radius(q) = sqrt(-λ / (π * imag(1 / q)))
        y0 = 100mm
        n1, n2 = 1.5, 1.7
        dl = doublet(n1, n2; y0)
        Msys = S(DOUBLET.r[3], n2, 1) * P(DOUBLET.l[2]) * S(DOUBLET.r[2], n1, n2) * P(DOUBLET.l[1]) *
               S(DOUBLET.r[1], 1, n1) * P(y0)
        ybk = y0 + sum(DOUBLET.l)
        for (kind, make) in (("GaussianBeamlet", () -> GaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0)),
            ("AstigmaticGaussianBeamlet", () -> AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, w0)))
            @testset "$kind through the doublet" begin
                g = make()
                solve_system!(System([dl]), g)
                comps = components(g)
                # every component beam resolves the same boundary sequence as the chief
                seq = [shape_sequence(c) for c in comps]
                @test length(comps) == (g isa GaussianBeamlet ? 3 : 9)
                @test all(s -> same_sequence(s, first(seq)), seq)
                @test length(first(seq)) == 4
                # front lens, then the cemented surface and the rear surface of the back lens
                @test same_sequence(first(seq), [BMO.shape(dl.front), BMO.shape(dl.back), BMO.shape(dl.back), nothing])
                # the chief follows the plain ray
                beam = solve([dl], beam_at([0, 0, 0], [0, 1, 0]; λ))
                @test flatten(g) ≈ flatten(beam) atol = TOL_POS
                # beam size behind the doublet against the ABCD law
                for z in (20mm, 60mm)
                    w = BMO.gauss_parameters(g, ybk + z)
                    w_ref = radius(prop(P(z) * Msys, q0))
                    @test all(x -> isapprox(x, w_ref; rtol = 1e-4), g isa GaussianBeamlet ? w[1:1] : w[1:2])
                end
                # and against the equivalent singlet for equal glass
                dl_eq = doublet(1.6, 1.6; y0)
                single = SphericalLens(DOUBLET.r[1], DOUBLET.r[3], sum(DOUBLET.l), DOUBLET.d, λ -> 1.6)
                translate3d!(single, [0, y0, 0])
                ga, gb = make(), make()
                solve_system!(System([dl_eq]), ga)
                solve_system!(System([single]), gb)
                for z in (20mm, 60mm)
                    @test all(BMO.gauss_parameters(ga, ybk + z)[1:2] .≈ BMO.gauss_parameters(gb, ybk + z)[1:2])
                end
            end
        end

        @testset "Triplet with the astigmatic beamlet: all rays resolve the same shapes" begin
            tl = triplet(1.5, 1.74, 1.62; y0)
            g = AstigmaticGaussianBeamlet([0.0, 0, 0], [0.0, 1, 0], λ, 0.2mm)
            solve_system!(System([tl]), g)
            seq = [shape_sequence(c) for c in components(g)]
            @test all(s -> same_sequence(s, first(seq)), seq)
            @test length(first(seq)) == 5
        end

        # plate: beam from y = -50 mm to the coating at y = 0, normal incidence
        for kind in ("GaussianBeamlet", "AstigmaticGaussianBeamlet")
            make(p) = kind == "GaussianBeamlet" ? GaussianBeamlet(p, [0.0, 1, 0], λ, w0) :
                      AstigmaticGaussianBeamlet(p, [0.0, 1, 0], λ, w0)
            wsize(b, z) = first(BMO.gauss_parameters(b, z))
            @testset "$kind with the plate beamsplitter" begin
                t, n, R = 3mm, 1.5, 0.3
                pbs = RectangularPlateBeamsplitter(36mm, 25mm, t, λ -> n; reflectance = R)
                g = make([0.0, -50mm, 0])
                solve_system!(System([pbs]), g)
                @test length(leaves(g)) == 2
                for b in all_beams(g)
                    seq = [shape_sequence(c) for c in components(b)]
                    @test all(s -> same_sequence(s, first(seq)), seq)
                end
                tr, rf = BMO.children(g)
                ref_beam = solve([pbs], beam_at([0, -50mm, 0], [0, 1, 0]; λ))
                @test flatten(g) ≈ flatten(ref_beam) atol = TOL_POS
                # q parameter: free space for the reflected beam, the plate shortens the effective distance by t (1 - 1/n)
                for z in (80mm, 120mm)
                    w_free = radius(q0 + z)
                    w_plate = radius(q0 + z - t * (1 - 1 / n))
                    @test wsize(rf, z) ≈ w_free rtol = 1e-4
                    @test wsize(tr, z) ≈ w_plate rtol = 1e-4
                end
                # field amplitudes behind the coating split as r : t
                # The field behind the coating splits as t : r. The field of an astigmatic beamlet is a physical field
                # (Poynting flux ~ n |E|²), such that it is smaller by sqrt(1 / n) in the glass of the plate. The field of a
                # GaussianBeamlet is not rescaled by the medium.
                f_medium = kind == "GaussianBeamlet" ? 1.0 : sqrt(1 / n)
                @test norm(BMO.electric_field(tr)) / norm(BMO.electric_field(rf)) ≈ f_medium * sqrt((1 - R) / R) rtol = 1e-9
            end

            @testset "$kind with the cube beamsplitter" begin
                L, n, R = 20mm, 1.5, 0.3
                cbs = CubeBeamsplitter(L, n; reflectance = R)
                g = make([0.0, -40mm, 0])
                solve_system!(System([cbs]), g)
                @test length(leaves(g)) == 2
                for b in all_beams(g)
                    seq = [shape_sequence(c) for c in components(b)]
                    @test all(s -> same_sequence(s, first(seq)), seq)
                end
                tr, rf = BMO.children(g)
                ref_beam = solve([cbs], beam_at([0, -40mm, 0], [0, 1, 0]; λ))
                @test flatten(g) ≈ flatten(ref_beam) atol = TOL_POS
                # both branches travel L in glass in total
                for z in (80mm, 120mm)
                    w_ref = radius(q0 + z - L * (1 - 1 / n))
                    @test wsize(tr, z) ≈ w_ref rtol = 1e-4
                    @test wsize(rf, z) ≈ w_ref rtol = 1e-4
                end
                @test norm(BMO.electric_field(tr)) / norm(BMO.electric_field(rf)) ≈ sqrt((1 - R) / R) rtol = 1e-9
            end
        end

        @testset "PolarizedRay: Fresnel chain through the cemented doublet and triplet" begin
            T(a, b) = 4a * b / (a + b)^2
            E0 = [0.6, 0.8, 0] # |E|² = 1; transverse for the direction +y is the x-z plane, use x and z
            E0 = [0.6, 0, 0.8]
            flux(ray) = BMO.refractive_index(ray) * sum(abs2, BMO.polarization(ray))
            for (obj, ns) in ((doublet(1.5, 1.72; y0 = 30mm), (1.5, 1.72)),
                (doublet(1.6, 1.6; y0 = 30mm), (1.6, 1.6)),
                (triplet(1.5, 1.74, 1.62; y0 = 30mm), (1.5, 1.74, 1.62)))
                beam = Beam([0.0, -20mm, 0], [0.0, 1, 0], λ0, E0)
                solve_system!(System([obj]), beam)
                chain = [1.0; collect(ns); 1.0]
                P_ref = prod(T(a, b) for (a, b) in zip(chain[1:(end - 1)], chain[2:end]))
                @test flux(last(rays(beam))) ≈ P_ref atol = 1e-12
                # no power change without an index step: same glass gives T = 1 at the cemented surface
            end
        end
    end
end


#=
The hits of a ray at one point, see `BMO.coincident_hits`: the solver collects all surfaces that a
ray hits within `get_coincidence_tolerance()` of the closest hit. This is checked directly on
`trace_all` and `trace_one`, for rays that start inside the first body.
=#

"""The hits of a ray from `p` along `d` in the `system`, traced like a solve does it."""
function hits_at(system, p, d; hint = nothing)
    ray = Ray(collect(Float64, p), collect(Float64, d), λ0)
    return BMO.with_bounding_spheres(system) do
        isnothing(hint) ? BMO.trace_all(system, ray) : BMO.trace_one(system, ray, hint)
    end
end

hit_shapes(hit) = [BMO.shape(h) for h in BMO.coincident_hits(hit)]
hit_lengths(hit) = [length(h) for h in BMO.coincident_hits(hit)]
holds(shapes, s) = any(x -> x === s, shapes)

@testset "Coincident hits are collected" begin
    tol = BMO.get_coincidence_tolerance()
    @test tol == 1e-9

    @testset "Doublet: the cemented faces of front and back" begin
        dl = doublet(1.5, 1.7)
        system = System([dl])
        # from within the front lens to the cemented surface at y = 4 mm
        hit = hits_at(system, [1mm, 2mm, 0], [0, 1, 0])
        shapes = hit_shapes(hit)
        @test length(shapes) == 2
        @test holds(shapes, BMO.shape(dl.front)) && holds(shapes, BMO.shape(dl.back))
        @test all(t -> abs(t - first(hit_lengths(hit))) ≤ tol, hit_lengths(hit))
        @test all(h -> BMO.object(h) === dl, BMO.coincident_hits(hit))
        # the outward normals of both faces are opposite
        n_front, n_back = (BMO.normal3d(h) for h in BMO.coincident_hits(hit) if BMO.shape(h) === BMO.shape(dl.front)),
            (BMO.normal3d(h) for h in BMO.coincident_hits(hit) if BMO.shape(h) === BMO.shape(dl.back))
        @test dot(first(n_front), first(n_back)) ≈ -1 atol = 1e-6
        # the front face of the doublet is a single hit
        @test length(hit_shapes(hits_at(system, [1mm, -10mm, 0], [0, 1, 0]))) == 1
    end

    @testset "Doublet with a hint: the hinted face first, the touching face added" begin
        dl = doublet(1.5, 1.7)
        system = System([dl])
        hit = hits_at(system, [1mm, 2mm, 0], [0, 1, 0]; hint = BMO.Hint(dl, BMO.shape(dl.back)))
        shapes = hit_shapes(hit)
        @test first(shapes) === BMO.shape(dl.back)
        @test length(shapes) == 2 && holds(shapes, BMO.shape(dl.front))
    end

    @testset "Triplet: two cemented surfaces, one at a time" begin
        tl = triplet(1.5, 1.7, 1.6)
        system = System([tl])
        v = triplet_vertices(0.0)
        hit = hits_at(system, [1mm, (v[2] + v[3]) / 2, 0], [0, 1, 0])
        shapes = hit_shapes(hit)
        @test length(shapes) == 2
        @test holds(shapes, BMO.shape(tl.middle)) && holds(shapes, BMO.shape(tl.back))
    end

    @testset "Cube beamsplitter: front, coating and back at the hypotenuse" begin
        cbs = CubeBeamsplitter(20mm, 1.5)
        system = System([cbs])
        hit = hits_at(system, [2mm, -8mm, 0], [0, 1, 0])
        shapes = hit_shapes(hit)
        @test length(shapes) == 3
        @test holds(shapes, BMO.shape(cbs.front)) && holds(shapes, BMO.shape(cbs.back)) &&
              holds(shapes, BMO.shape(cbs.coating))
    end

    @testset "Separate touching objects, also within a group" begin
        A = RectangularCompensatorPlate(40mm, 40mm, 5mm, λ -> 1.7)
        B = RectangularCompensatorPlate(40mm, 40mm, 3mm, λ -> 1.4)
        translate3d!(B, [0, 5mm, 0])
        d = normalize([0.3, 1, 0.1])
        for system in (System([A, B]), System([B, A]), System(ObjectGroup([A, B])), StaticSystem([A, B]))
            hit = hits_at(system, [0, 2mm, 0], d)
            objs = [BMO.object(h) for h in BMO.coincident_hits(hit)]
            @test length(objs) == 2
            @test any(o -> o === A, objs) && any(o -> o === B, objs)
        end
    end

    @testset "Gaps: below the tolerance a contact, above two hits" begin
        for (gap, n_hits) in ((0.3e-9, 2), (5e-9, 1), (1e-6, 1))
            A = RectangularCompensatorPlate(40mm, 40mm, 5mm, λ -> 1.7)
            B = RectangularCompensatorPlate(40mm, 40mm, 3mm, λ -> 1.4)
            translate3d!(B, [0, 5mm + gap, 0])
            hit = hits_at(System([A, B]), [0, 2mm, 0], [0, 1, 0])
            @test length(hit_shapes(hit)) == n_hits
            @test BMO.object(hit) === A
        end
    end

    @testset "Skipping parts beyond the closest hit gives the same hits" begin
        # lenses and plates scattered in space, rays from the origin in many directions: the pruned
        # trace equals the closest hits of all objects without bounding spheres
        objs = Any[]
        for k in 1:12
            o = isodd(k) ? SphericalLens(30mm, -30mm, 5mm, 15mm, λ -> 1.5) :
                RectangularCompensatorPlate(10mm, 10mm, 2mm, λ -> 1.6)
            φ = 2π * k / 12
            translate3d!(o, [60mm * cos(φ), 60mm * sin(φ), 10mm * sin(3φ)])
            push!(objs, o)
        end
        # a touching pair among them
        A = RectangularCompensatorPlate(10mm, 10mm, 2mm, λ -> 1.6)
        B = RectangularCompensatorPlate(10mm, 10mm, 2mm, λ -> 1.4)
        translate3d!(A, [0, 30mm, 0])
        translate3d!(B, [0, 32mm, 0])
        append!(objs, [A, B])
        system = System(objs)
        for k in 0:359
            φ = deg2rad(k + 0.5)
            for p in ([0.0, 0, 0], [0.0, 31mm, 0])
                ray = Ray(p, [cos(φ), sin(φ), 0.02], λ0)
                pruned = BMO.with_bounding_spheres(() -> BMO.trace_all(system, ray), system)
                plain = BMO.intersect3d(objs, ray)
                @test isnothing(pruned) == isnothing(plain)
                isnothing(plain) && continue
                @test length(pruned) == length(plain)
                @test Set(objectid.(hit_shapes(pruned))) == Set(objectid.(hit_shapes(plain)))
            end
        end
    end
end

end # module
