module TestStorageObjects

using BeamletOptics
using Test
using LinearAlgebra: normalize, norm

const BMO = BeamletOptics

const mm = 1e-3
const inch = 25.4mm
const NBK7 = SellmeierEquation(1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653)

tmpfile() = joinpath(mktempdir(), "objects.bmo")

"Saves `objs` as one system, loads it and returns the loaded objects"
function roundtrip(objs...; names = IdDict{Any, String}())
    path = tmpfile()
    save_setup(path; systems = [System(collect(BMO.AbstractObject, objs))], names)
    return Tuple(load_setup(path).systems[1].objects)
end

"Beams of the traced tree in depth-first order"
function beam_tree(beam, out = Any[])
    push!(out, beam)
    foreach(c -> beam_tree(c, out), beam.children)
    return out
end

function trace(obj, ray)
    beam = Beam(deepcopy(ray))
    solve_system!(System([obj]), beam)
    return beam_tree(beam)
end

"Compares the traced beam trees: same number of rays, positions within 1e-10 and directions within 1e-12"
function same_trace(t1, t2)
    length(t1) == length(t2) || return false
    for (b1, b2) in zip(t1, t2)
        r1, r2 = BMO.rays(b1), BMO.rays(b2)
        length(r1) == length(r2) || return false
        for (a, b) in zip(r1, r2)
            norm(position(a) - position(b)) < 1e-10 || return false
            norm(direction(a) - direction(b)) < 1e-12 || return false
            if a isa PolarizedRay
                norm(BMO.polarization(a) - BMO.polarization(b)) < 1e-10 || return false
            end
        end
    end
    return true
end

"Moves `obj` to a non-trivial pose"
function place!(obj)
    translate3d!(obj, [3mm, 50mm, -2mm])
    zrotate3d!(obj, 0.2)
    xrotate3d!(obj, 0.1)
    return obj
end

"Test ray that starts in front of `obj` and runs roughly along +y through its center"
function test_ray(obj; polarized = false)
    dir = normalize([0.02, 1.0, 0.01])
    pos = position(obj) - 0.2 * dir + [0.5mm, 0, 0.5mm]
    E0 = [1.0, 0, 1.0] - ([1.0, 0, 1.0]' * dir) * dir
    return polarized ? PolarizedRay(pos, dir, 1e-6, E0) : Ray(pos, dir, 1e-6)
end

@testset "Storage of objects and materials" begin
    @testset "ConstantRefractiveIndex" begin
        c = ConstantRefractiveIndex(1.5)
        @test c(532e-9) == 1.5
        @test c(1) == 1.5
        @test c isa BMO.RefractiveIndex
        @test isnothing(BMO.test_refractive_index_function(c))
        @test SphericalLens(50mm, -50mm, 5mm).n isa ConstantRefractiveIndex
        @test SphericalLens(50mm, -50mm, 5mm, 1inch, 1.6).n == ConstantRefractiveIndex(1.6)
        @test ThinLens(50mm, -50mm, 1inch, 1.6).n == ConstantRefractiveIndex(1.6)
        @test RectangularCompensatorPlate(10mm, 10mm, 2mm, 1.6).n == ConstantRefractiveIndex(1.6)
    end

    @testset "Materials" begin
        ctx = BMO.StorageContext(BMO.bmo_version())
        rt(n) = BMO.decode_material(BMO.encode_material(n, ctx), ctx)
        @test rt(ConstantRefractiveIndex(1.4567)) === ConstantRefractiveIndex(1.4567)
        @test BMO.encode_material(ConstantRefractiveIndex(1.5), ctx) == Dict("type" => "Constant", "n" => 1.5)
        @test rt(NBK7) === NBK7
        @test BMO.encode_material(NBK7, ctx)["type"] == "Sellmeier"
        dri = DiscreteRefractiveIndex([1064e-9, 532e-9, 633e-9], [1.50, 1.52, 1.51])
        loaded = rt(dri)
        @test loaded isa DiscreteRefractiveIndex{Float64}
        @test loaded.data == dri.data
        @test BMO.encode_material(dri, ctx)["lambda"]["data"] == [532e-9, 633e-9, 1064e-9]

        # Registered closures are stored by name
        f = λ -> 1.5 + 1e-15 / λ^2
        @test_throws r"register_material!" BMO.encode_material(f, ctx)
        @test_throws ArgumentError BMO.encode_material(f, ctx)
        register_material!("test-glass", f)
        d = BMO.encode_material(f, ctx)
        @test d == Dict("type" => "Named", "name" => "test-glass")
        @test rt(f) === f
        @test_throws r"material \"unknown-glass\" is not registered" BMO.decode_material(
            Dict{String, Any}("type" => "Named", "name" => "unknown-glass"), ctx)
        @test_throws ArgumentError register_material!("", f)
        @test_throws ArgumentError register_material!("not-callable", λ -> "glass")
        delete!(BMO.MATERIALS, "test-glass")
    end

    @testset "Unregistered closure" begin
        lens = SphericalLens(50mm, -50mm, 5mm, 1inch, λ -> 1.5)
        path = tmpfile()
        err = try
            save_setup(path; systems = [System([lens])], names = IdDict{Any, String}(lens => "my-lens"))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("register_material!", err.msg)
        @test occursin("Lens", err.msg)
        # TODO: the object id is only known in Core.jl (`_collect_object!`), which must add it to the message
        @test occursin("my-lens", err.msg)
        @test !isfile(path)

        # A registered closure round-trips
        n = λ -> 1.55
        register_material!("closure-glass", n)
        lens = SphericalLens(50mm, -50mm, 5mm, 1inch, n)
        loaded, = roundtrip(lens)
        @test loaded.n === n
        delete!(BMO.MATERIALS, "closure-glass")
    end

    # Every exported component constructor, name => (object, polarized test ray)
    components = [
        "SquarePlanoMirror2D" => SquarePlanoMirror2D(20mm),
        "RectangularPlanoMirror" => RectangularPlanoMirror(20mm, 10mm, 5mm),
        "SquarePlanoMirror" => SquarePlanoMirror(20mm, 5mm),
        "RoundPlanoMirror" => RoundPlanoMirror(1inch, 5mm),
        "SphericalMirror" => SphericalMirror(100mm, 10mm, 1inch),
        "RightAnglePrismMirror" => RightAnglePrismMirror(20mm, 20mm),
        "OffAxisParabolicMirror" => OffAxisParabolicMirror(50mm, 1inch; angle = 60),
        "ParabolicMirror" => ParabolicMirror(50mm, 1inch),
        "SphericalLens" => SphericalLens(50mm, -80mm, 5mm, 1inch, NBK7),
        "SphericalLens (default n)" => SphericalLens(50mm, Inf, 5mm),
        "ThinLens" => ThinLens(50mm, -50mm, 1inch, 1.5),
        "SphericalDoubletLens" => SphericalDoubletLens(60mm, -40mm, -200mm, 6mm, 3mm, 1inch, NBK7,
            DiscreteRefractiveIndex([1e-6], [1.62])),
        "SphericalTripletLens" => SphericalTripletLens(60mm, -40mm, 40mm, -60mm, 5mm, 2mm, 5mm, 1inch,
            NBK7, ConstantRefractiveIndex(1.7), NBK7),
        "Lens (aspheric)" => Lens(EvenAsphericalSurface(20mm, 1inch, -1.0, [0, 2.1e2, 1.7e1]),
            SphericalSurface(-100mm, 1inch), 8mm, NBK7),
        "Lens (flat)" => Lens(CircularFlatSurface(1inch), 3mm, ConstantRefractiveIndex(1.5)),
        "Lens (rectangular flat)" => Lens(RectangularFlatSurface(20mm),
            RectangularFlatSurface(20mm), 3mm, ConstantRefractiveIndex(1.5)),
        "Lens (cylindrical)" => Lens(CylindricalSurface(20mm, 10mm, 20mm), 5mm, ConstantRefractiveIndex(1.517)),
        "Lens (acylindrical)" => Lens(AcylindricalSurface(15mm, 25mm, 50mm, -1.0, [0, 1.19e1, -2.9e3]),
            7.5mm, ConstantRefractiveIndex(1.777)),
        "RightAnglePrism" => RightAnglePrism(20mm, 20mm, NBK7),
        "Detector" => Detector(30mm),
        "Detector (no stop)" => Detector(30mm, false),
        "ThinBeamsplitter" => ThinBeamsplitter(20mm, 20mm; reflectance = 0.3),
        "RoundThinBeamsplitter" => RoundThinBeamsplitter(1inch; reflectance = 0.7),
        "RectangularPlateBeamsplitter" => RectangularPlateBeamsplitter(30mm, 30mm, 3mm, NBK7; reflectance = 0.4),
        "RoundPlateBeamsplitter" => RoundPlateBeamsplitter(1inch, 3mm, NBK7),
        "CubeBeamsplitter" => CubeBeamsplitter(20mm, NBK7; reflectance = 0.6),
        "RectangularCompensatorPlate" => RectangularCompensatorPlate(30mm, 30mm, 3mm, 1.5),
        "PolarizationFilter" => PolarizationFilter(30mm),
        "RoundPolarizationFilter" => RoundPolarizationFilter(1inch),
        "RoundLinearPolarizer" => RoundLinearPolarizer(1inch, 1.6mm, 1.6mm, NBK7),
        "NonInteractableObject" => NonInteractableObject(BMO.CuboidMesh(10mm, 10mm, 10mm)),
        "IntersectableObject" => IntersectableObject(BMO.CuboidMesh(10mm, 10mm, 10mm)),
        "Retroreflector" => Retroreflector(0.02),
    ]

    @testset "Component round trip: $name" for (name, obj) in components
        place!(obj)
        loaded, = roundtrip(obj)
        @test typeof(loaded) == typeof(obj)
        @test position(loaded) ≈ position(obj) atol = 1e-14
        @test orientation(loaded) ≈ orientation(obj) atol = 1e-14
        polarizer = obj isa Union{PolarizationFilter, LinearPolarizer}
        for polarized in unique((false, polarizer))
            ray = test_ray(obj; polarized)
            t1, t2 = trace(obj, ray), trace(loaded, ray)
            @test same_trace(t1, t2)
            # The test ray hits the object
            obj isa NonInteractableObject || @test !isnothing(BMO.intersection(first(BMO.rays(t1[1]))))
        end
    end

    @testset "Component fields" begin
        bs = ThinBeamsplitter(20mm, 20mm; reflectance = 0.3)
        loaded, = roundtrip(bs)
        @test loaded.reflectance == bs.reflectance
        @test loaded.transmittance == bs.transmittance

        pf = place!(PolarizationFilter(30mm; cutoff_strength = 1e-6))
        loaded, = roundtrip(pf)
        @test loaded.JMat.data == pf.JMat.data
        @test loaded.cutoff == pf.cutoff

        # Hits are not stored, the loaded detector has its own lock
        det = place!(Detector(30mm))
        solve_system!(System([det]), Beam(test_ray(det)))
        @test !isnothing(BMO.hits(det))
        loaded, = roundtrip(det)
        @test isnothing(BMO.hits(loaded))
        @test loaded.stop == det.stop
        @test loaded.lock !== det.lock

        # Sub-objects are stored inline, the loaded system has only the top-level object
        cbs = CubeBeamsplitter(20mm, NBK7)
        loaded, = roundtrip(cbs)
        @test loaded.front.n === NBK7
        @test loaded.coating.reflectance == cbs.coating.reflectance

        # An object that occurs in several systems is stored once
        m = SquarePlanoMirror(20mm, 5mm)
        path = tmpfile()
        save_setup(path; systems = [System([m]), System([m])])
        setup = load_setup(path)
        @test setup.systems[1].objects[1] === setup.systems[2].objects[1]
    end
end

end # module
