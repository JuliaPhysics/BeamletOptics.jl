module TestSurfaceStacks

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 550e-9

flux(r) = BMO.refractive_index(r) * norm(BMO.polarization(r))^2
polarized(E; p = [0, 0, 0.0], d = [0, 1.0, 0]) = Beam(PolarizedRay(p, d, λ0, complex.(E)))

# A free-standing film at y = 50 mm, its front towards the source (-y)
function film(model)
    layer = CoatingLayer(BMO.RectangularFlatMesh(20mm, 20mm), model)
    zrotate3d!(layer, π)
    translate3d!(layer, [0, 50mm, 0])
    return layer
end

function solve_film(model, E; kwargs...)
    beam = polarized(E)
    solve_system!(System([film(model)]), beam; kwargs...)
    return beam
end

# The power that leaves towards +y and -y, summed over the leaves of the beam tree
function exit_powers(beam, y)
    forward, backward = 0.0, 0.0
    for leaf in Leaves(beam)
        r = last(rays(leaf))
        BMO.direction(r)[2] > 0 && BMO.position(r)[2] > y - 1e-6 && (forward += flux(r))
        BMO.direction(r)[2] < 0 && (backward += flux(r))
    end
    return forward, backward
end

@testset "Retarders" begin
    # A quarter-wave film at 45° turns linear into circular polarization
    E = BMO.polarization(rays(solve_film(Retarder(π / 2, [1, 0, 1]), [1, 0, 0]))[2])
    @test abs(E[1]) ≈ abs(E[3]) ≈ 1 / √2
    @test abs(angle(E[3] / E[1])) ≈ π / 2
    # A half-wave film at 22.5° turns the polarization by 45°
    E = BMO.polarization(rays(solve_film(Retarder(π, [cos(π / 8), 0, sin(π / 8)]), [1, 0, 0]))[2])
    @test abs.(E) ≈ [1, 0, 1] / √2 atol = 1e-12
    # A zero-order plate is a half-wave plate at half its design wavelength
    qwp = Retarder(π / 2, [1, 0, 0]; λ0 = 1100e-9)
    @test BMO.retardance(qwp, λ0) ≈ π
    @test BMO.retardance(Retarder(π / 2, [1, 0, 0]), λ0) == π / 2
end

@testset "Reflective polarizers" begin
    # Malus's law in transmission, the rest reflected, at any orientation of the film
    for φ in (0, π / 6, π / 3, π / 2)
        beam = solve_film(ReflectivePolarizer([cos(φ), 0, sin(φ)]), [1, 0, 0])
        forward, backward = exit_powers(beam, 50mm)
        @test forward ≈ cos(φ)^2 atol = 1e-12
        @test backward ≈ sin(φ)^2 atol = 1e-12
    end
    # The transmitted field lies along the axis
    beam = solve_film(ReflectivePolarizer([0, 0, 1]), [1, 0, 1] / √2)
    @test abs.(BMO.polarization(first(rays(children(beam)[1])))) ≈ [0, 0, 1 / √2] atol = 1e-12
end

@testset "Stacks" begin
    # Two lossless 50:50 reflectors without a gap: R = 4R₁/(1 + R₁)², Born & Wolf §7.6
    bs = BeamsplitterCoating()
    forward, backward = exit_powers(solve_film(SurfaceStack(bs, bs), [1, 0, 0]), 50mm)
    @test backward ≈ 4 * 0.5 / 1.5^2
    @test forward + backward ≈ 1
    # The order of the elements matters: a polarizer that the ray meets first blocks x, behind a
    # quarter-wave film half of the circular light passes
    qwp, pol = Retarder(π / 2, [1, 0, 1]), ReflectivePolarizer([0, 0, 1])
    forward, _ = exit_powers(solve_film(SurfaceStack(qwp, pol), [1, 0, 0]), 50mm)
    @test forward ≈ 0 atol = 1e-12
    forward, _ = exit_powers(solve_film(SurfaceStack(pol, qwp), [1, 0, 0]), 50mm)
    @test forward ≈ 0.5
    # A stack of one element is the element
    forward, backward = exit_powers(solve_film(SurfaceStack(bs), [1, 0, 0]), 50mm)
    @test (forward, backward) == exit_powers(solve_film(bs, [1, 0, 0]), 50mm)
end

# A flat pancake fold: a 50:50 mirror on the front face of a plate, a quarter-wave film under a
# reflective polarizer on its back face. Circular light passes once, is reflected by the polarizer,
# reflected by the mirror with its handedness flipped and passes the polarizer: 1/4 of the input.
function pancake(Γ = π / 2; E = [1, 0, im] / √2)
    plate = RectangularCompensatorPlate(30mm, 30mm, 5mm, ConstantRefractiveIndex(1.5))
    surface_model!(plate, FlatFace([0, -1, 0]), BeamsplitterCoating())
    surface_model!(plate, FlatFace([0, 1, 0]),
        SurfaceStack(Retarder(Γ, [1, 0, -1]), ReflectivePolarizer([1, 0, 0])))
    zrotate3d!(plate, deg2rad(20))
    translate3d!(plate, [3mm, 50mm, 0])
    d = orientation(plate)[:, 2]
    beam = Beam(PolarizedRay(position(plate) - 40mm * d, d, λ0,
        orientation(plate) * complex.(E)))
    solve_system!(System([plate]), beam; depth_max = 30)
    forward, backward = 0.0, 0.0
    for leaf in Leaves(beam)
        r = last(rays(leaf))
        s = dot(BMO.direction(r), d)
        s > 0 && dot(BMO.position(r) - position(plate), d) > 0 && (forward += flux(r))
        s < 0 && (backward += flux(r))
    end
    return forward, backward, beam
end

@testset "A pancake fold" begin
    forward, backward, _ = pancake()
    @test forward ≈ 0.25
    @test backward ≈ 0.75
    # The other handedness passes directly: half of it
    forward, backward, _ = pancake(; E = [1, 0, -im] / √2)
    @test forward ≈ 0.5
    @test forward + backward ≈ 1
    # A retardance error δ leaks the direct path, the ghost: of the half that enters, sin²(δ/2) passes
    # the polarizer at once
    for δ in (0.05, 0.2)
        forward, backward, beam = pancake(π / 2 + δ)
        @test forward + backward ≈ 1
        ghost = first(rays(children(children(beam)[1])[1]))
        @test flux(ghost) ≈ 0.5 * sin(δ / 2)^2
    end
end

end # module
