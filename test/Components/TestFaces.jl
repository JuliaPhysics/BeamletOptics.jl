module TestFaces

using BeamletOptics
using Test
using LinearAlgebra
using AbstractTrees

const BMO = BeamletOptics
const mm = 1e-3
const λ0 = 550e-9
const n_glass = 1.5
const glass = ConstantRefractiveIndex(n_glass)

# Power flux of a polarized ray, relative to an input field of unit norm in air
flux(r) = BMO.refractive_index(r) * norm(BMO.polarization(r))^2
fresnel_loss = ((n_glass - 1) / (n_glass + 1))^2

polarized(p, d, E = [1.0 + 0im, 0, 0]) = Beam(PolarizedRay(p, d, λ0, E))

function lens_at(p)
    lens = SphericalLens(50mm, -50mm, 5mm, 25mm, glass)
    translate3d!(lens, p)
    return lens
end

@testset "Coatings on the faces of a lens" begin
    lens = lens_at([0, 50mm, 0])
    system = System([lens])

    beam = polarized([0, 0, 0], [0, 1, 0])
    solve_system!(system, beam)
    @test flux(rays(beam)[2]) ≈ 1 - fresnel_loss
    @test flux(rays(beam)[3]) ≈ (1 - fresnel_loss)^2

    # The front face only
    surface_model!(lens, FrontFace(), AntireflectionCoating())
    @test surface_model(lens, FrontFace()) isa BeamsplitterCoating
    @test isnothing(surface_model(lens, BackFace()))
    beam = polarized([0, 0, 0], [0, 1, 0])
    solve_system!(system, beam)
    @test flux(rays(beam)[2]) ≈ 1
    @test flux(rays(beam)[3]) ≈ 1 - fresnel_loss

    # Both faces, and the coating moves with the lens
    surface_model!(lens, BackFace(), AntireflectionCoating())
    zrotate3d!(lens, deg2rad(30))
    translate3d!(lens, [5mm, 0, 2mm])
    d, E = orientation(lens)[:, 2], complex.(orientation(lens)[:, 1])
    beam = polarized(position(lens) - 40mm * d, d, E)
    solve_system!(system, beam)
    @test flux(rays(beam)[3]) ≈ 1

    # Removing a coating
    surface_model!(lens, FrontFace(), nothing)
    @test isnothing(surface_model(lens, FrontFace()))
    beam = polarized(position(lens) - 40mm * d, d, E)
    solve_system!(system, beam)
    @test flux(rays(beam)[3]) ≈ 1 - fresnel_loss

    # The edge is a face of its own
    @test isnothing(surface_model(lens, EdgeFace()))

    # A part whose faces carry no models
    @test_throws ArgumentError surface_model!(SquarePlanoMirror2D(10mm), AllFaces(), MirrorCoating())
end

@testset "A mask on a face" begin
    lens = lens_at([0, 50mm, 0])
    # A mirror within r < 5 mm of the axis of the lens, a mask along its local y-axis
    surface_model!(lens, MaskedFace(FrontFace(), BMO.CylinderSDF(5mm, 20mm)), MirrorCoating())
    system = System([lens])
    beam = Beam([2mm, 0, 0], [0, 1, 0], λ0)
    solve_system!(system, beam)
    @test BMO.direction(rays(beam)[2])[2] < 0
    beam = Beam([8mm, 0, 0], [0, 1, 0], λ0)
    solve_system!(system, beam)
    @test BMO.direction(rays(beam)[2])[2] > 0
end

@testset "Flat faces of a plate" begin
    plate = RectangularCompensatorPlate(40mm, 40mm, 5mm, glass)
    translate3d!(plate, [0, 50mm, 0])
    surface_model!(plate, FlatFace([0, -1, 0]), MirrorCoating())
    beam = Beam([0, 0, 0], [0, 1, 0], λ0)
    solve_system!(System([plate]), beam)
    @test BMO.direction(last(rays(beam))) ≈ [0, -1, 0]
    # The opposite face is not coated
    beam = Beam([0, 100mm, 0], [0, -1, 0], λ0)
    solve_system!(System([plate]), beam)
    @test BMO.direction(rays(beam)[2]) ≈ [0, -1, 0]
end

@testset "A coating at a cemented face lies on its own body" begin
    # A thin silver film on the back face of plate A, cemented to plate B. Absorbing films reflect
    # differently from both sides, hence the side of the substrate matters.
    silver = AbsorbingRefractiveIndex(ConstantRefractiveIndex(0.135), 3.985)
    film = ThinFilmCoating(ThinFilmLayer(silver, 20e-9))
    nA, nB, tA = 1.5, 1.7, 5mm
    A = RectangularCompensatorPlate(40mm, 40mm, tA, ConstantRefractiveIndex(nA))
    B = RectangularCompensatorPlate(40mm, 40mm, 5mm, ConstantRefractiveIndex(nB))
    translate3d!(B, [0, tA, 0])
    surface_model!(A, FlatFace([0, 1, 0]), film)
    for (beam, n_in, n_out, layers) in (
            (polarized([0, -10mm, 0], [0, 1, 0]), nA, nB, film.layers),
            (polarized([0, 20mm, 0], [0, -1, 0]), nB, nA, reverse(film.layers)))
        solve_system!(System([A, B]), beam; depth_max = 2)
        _, (R, _, _, _) = BMO.thin_film_coefficients(layers, n_in, n_out, 1.0, λ0)
        incoming = last(rays(beam))
        reflected = first(rays(children(beam)[2]))
        @test flux(reflected) / flux(incoming) ≈ R
    end
    _, (R_A, _, _, _) = BMO.thin_film_coefficients(film.layers, nA, nB, 1.0, λ0)
    _, (R_B, _, _, _) = BMO.thin_film_coefficients(reverse(film.layers), nB, nA, 1.0, λ0)
    @test !isapprox(R_A, R_B; rtol = 1e-3)
end

# Two plates cemented at y = 5 mm, A below and B above
function cemented_plates()
    A = RectangularCompensatorPlate(40mm, 40mm, 5mm, glass)
    B = RectangularCompensatorPlate(40mm, 40mm, 5mm, glass)
    translate3d!(B, [0, 5mm, 0])
    return A, B
end

@testset "Two coated faces of a cemented boundary act together" begin
    # A mirror on A and an anti-reflection coating on B: both lie in the boundary, the light is
    # reflected from either side
    A, B = cemented_plates()
    surface_model!(A, FlatFace([0, 1, 0]), MirrorCoating())
    surface_model!(B, FlatFace([0, -1, 0]), AntireflectionCoating())
    for (p, d) in (([0, -10mm, 0], [0, 1, 0]), ([0, 20mm, 0], [0, -1, 0]))
        beam = Beam(p, d, λ0)
        solve_system!(System([A, B]), beam)
        @test BMO.direction(last(rays(beam))) ≈ -d
    end
    # Two thin films act as one film of all their layers: a silver film split between both faces
    silver = AbsorbingRefractiveIndex(ConstantRefractiveIndex(0.135), 3.985)
    half, whole = ThinFilmLayer(silver, 10e-9), ThinFilmLayer(silver, 20e-9)
    function reflected_power(models)
        A, B = cemented_plates()
        for (plate, normal, model) in zip((A, B), ([0, 1, 0], [0, -1, 0]), models)
            isnothing(model) || surface_model!(plate, FlatFace(normal), model)
        end
        beam = polarized([0, -10mm, 0], [0, 1, 0])
        solve_system!(System([A, B]), beam; depth_max = 2)
        return flux(first(rays(children(beam)[2]))) / flux(last(rays(beam)))
    end
    @test reflected_power((ThinFilmCoating(half), ThinFilmCoating(half))) ≈
          reflected_power((ThinFilmCoating(whole), nothing))
    # A retarder on the entered body acts in the frame of that body: a half-wave film on B, turned
    # with B, turns the light as the same film alone does
    function turned_field(coat_A)
        A, B = cemented_plates()
        coat_A && surface_model!(A, FlatFace([0, 1, 0]), AntireflectionCoating())
        surface_model!(B, FlatFace([0, -1, 0]), Retarder(π, [1, 0, 0]))
        rotate3d!(B, [0, 1, 0], deg2rad(22.5))
        beam = polarized([0, -10mm, 0], [0, 1, 0])
        solve_system!(System([A, B]), beam)
        return abs.(BMO.polarization(last(rays(beam))))
    end
    @test turned_field(true) ≈ turned_field(false)
    @test turned_field(true)[3] > 0.1
    # A pair that cannot be combined warns once, the coating of the body that the ray leaves acts
    A, B = cemented_plates()
    surface_model!(A, FlatFace([0, 1, 0]), MirrorCoating())
    surface_model!(B, FlatFace([0, -1, 0]), ThinFilmCoating(half))
    beam = Beam([0, -10mm, 0], [0, 1, 0], λ0)
    @test_logs (:warn, r"cannot be combined") solve_system!(System([A, B]), beam)
    @test BMO.direction(last(rays(beam))) ≈ [0, -1, 0]
end

@testset "The surface model of a system" begin
    lens = lens_at([0, 50mm, 0])
    beam = polarized([0, 0, 0], [0, 1, 0])
    solve_system!(System([lens]; surface_model = FresnelGhosts()), beam; depth_max = 2)
    @test length(children(beam)) == 2
    @test BMO.power_fraction(children(beam)[2], beam) ≈ fresnel_loss
    # A coated face takes precedence over the model of the system
    surface_model!(lens, FrontFace(), AntireflectionCoating())
    beam = polarized([0, 0, 0], [0, 1, 0])
    solve_system!(StaticSystem([lens]; surface_model = FresnelGhosts()), beam; depth_max = 2)
    @test flux(rays(beam)[2]) ≈ 1
    @test_throws ArgumentError System([lens]; surface_model = 1.0)
end

# A dummy face, e.g. of a shape that knows its faces natively: one quadrant of a plate
struct QuadrantFace <: BMO.AbstractFace end
BMO.isonface(::QuadrantFace, ::BMO.AbstractObject, ::BMO.Intersection, point, normal) =
    point[1] > 0 && point[3] > 0

@testset "A face needs only isonface" begin
    plate = RectangularCompensatorPlate(40mm, 40mm, 5mm, glass)
    translate3d!(plate, [0, 50mm, 0])
    surface_model!(plate, QuadrantFace(), MirrorCoating())
    system = System([plate])
    for (x, z, back) in ((5mm, 5mm, true), (-5mm, 5mm, false), (5mm, -5mm, false))
        beam = Beam([x, 0, z], [0, 1, 0], λ0)
        solve_system!(system, beam)
        @test (BMO.direction(rays(beam)[2])[2] < 0) == back
    end
end

# A dummy shape that knows its faces natively, like a boundary representation: a cuboid mesh whose
# hits carry the face that was hit
struct NativeFace <: BMO.AbstractFace
    id::Int
end

struct NativeFacesShape{M <: BMO.AbstractShape{Float64}} <: BMO.AbstractShape{Float64}
    mesh::M
end

BMO.position(s::NativeFacesShape) = position(s.mesh)
BMO.orientation(s::NativeFacesShape) = orientation(s.mesh)
BMO.translate3d!(::BMO.Movable, s::NativeFacesShape, offset) = translate3d!(s.mesh, offset)
BMO.rotate3d!(::BMO.Movable, s::NativeFacesShape, R::AbstractMatrix) = BMO.rotate3d!(s.mesh, R)

function BMO.intersect3d(s::NativeFacesShape, ray::BMO.AbstractRay)
    hit = BMO.intersect3d(s.mesh, ray)
    isnothing(hit) && return nothing
    n = BMO.normal3d(hit)
    # Faces numbered by their outward normal in the local frame: 1 front (-y), 2 back (+y), 3 others
    m = transpose(orientation(s)) * n
    id = m[2] < -0.5 ? 1 : m[2] > 0.5 ? 2 : 3
    return BMO.face!(BMO.Intersection(hit.t, n, s), NativeFace(id))
end

@testset "Native faces of a shape" begin
    mesh = BMO.CuboidMesh(40mm, 5mm, 40mm)
    translate3d!(mesh, [-20mm, 0, -20mm])
    BMO.set_new_origin3d!(mesh)
    body = Prism(NativeFacesShape(mesh), glass)
    translate3d!(body, [0, 50mm, 0])
    surface_model!(body, NativeFace(1), MirrorCoating())
    system = System([body])
    beam = Beam([0, 0, 0], [0, 1, 0], λ0)
    solve_system!(system, beam)
    @test BMO.direction(rays(beam)[2]) ≈ [0, -1, 0]
    # The back face is not coated; geometric faces and masks work on native faces as well
    beam = Beam([0, 100mm, 0], [0, -1, 0], λ0)
    solve_system!(system, beam)
    @test BMO.direction(rays(beam)[2]) ≈ [0, -1, 0]
    surface_model!(body, MaskedFace(NativeFace(2), BMO.CylinderSDF(5mm, 20mm)), MirrorCoating())
    beam = Beam([2mm, 100mm, 0], [0, -1, 0], λ0)
    solve_system!(system, beam)
    @test BMO.direction(rays(beam)[2]) ≈ [0, 1, 0]
end

end # module
