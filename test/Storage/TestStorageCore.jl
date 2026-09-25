module TestStorageCore

using BeamletOptics
using Test
using LinearAlgebra: I

const BMO = BeamletOptics

# Test-only types that use the extension hook
mutable struct StoreObj <: BMO.AbstractObject{Float64}
    pos::BMO.Point3{Float64}
end
BMO.to_storage(o::StoreObj, ctx) = Dict{String, Any}("pos" => BMO.encode_vec3(o.pos))
BMO.from_storage(::Type{StoreObj}, d, ctx) = StoreObj(BMO.decode_vec3(d["pos"]))

mutable struct StoreArrayObj <: BMO.AbstractObject{Float64}
    data::Matrix{Float64}
end
BMO.to_storage(o::StoreArrayObj, ctx) = Dict{String, Any}("data" => BMO.encode_array(o.data, ctx; stem = "data"))
BMO.from_storage(::Type{StoreArrayObj}, d, ctx) = StoreArrayObj(BMO.decode_array(d["data"], ctx))

mutable struct StoreSrc
    λ::Float64
end
BMO.to_storage(s::StoreSrc, ctx) = Dict{String, Any}("lambda" => s.λ)
BMO.from_storage(::Type{StoreSrc}, d, ctx) = StoreSrc(d["lambda"])

mutable struct Unregistered <: BMO.AbstractObject{Float64} end

BMO.register_storage_type!(StoreObj, "StoreObj")
BMO.register_storage_type!(StoreArrayObj, "StoreArrayObj")
BMO.register_storage_type!(StoreSrc, "StoreSrc")

tmpfile() = joinpath(mktempdir(), "setup.bmo")

@testset "Storage core" begin
    @testset "Version check" begin
        check(w, l) = BMO.check_version(VersionNumber(w); loader = VersionNumber(l))
        @test isnothing(check("0.14.0", "0.14.3"))
        @test_throws r"update BeamletOptics to ≥ 0.14.3" check("0.14.3", "0.14.0")
        @test_throws r"BeamletOptics 0.14.x release.*BeamletOptics@0.14" check("0.14.3", "0.15.0")
        @test isnothing(check("0.14.0", "0.14.0-DEV"))
        @test isnothing(check("1.0.0", "1.4.2"))
        @test_throws r"BeamletOptics 1.x release" check("1.9.0", "2.0.0")
        @test BMO.bmo_version() == pkgversion(BeamletOptics)
    end

    @testset "Array and value encoding" begin
        ctx = BMO.StorageContext(v"0.0.1")
        for A in (rand(3, 4), rand(1:10, 5), rand(ComplexF64, 2, 2), rand(3, 1000), rand(ComplexF64, 600))
            d = BMO.encode_array(A, ctx)
            @test haskey(d, "asset") == (length(A) > BMO.INLINE_ARRAY_LIMIT)
            B = BMO.decode_array(d, ctx)
            @test B == A
            @test typeof(B) == typeof(collect(A))
        end
        @test_throws ArgumentError BMO.encode_array(rand(Float32, 3), ctx)
        M = BMO.SMatrix{3, 3}(rand(3, 3))
        @test BMO.decode_mat3(BMO.encode_mat3(M)) == M
        @test BMO.encode_mat3(M)[2][3] == M[2, 3]
        @test BMO.decode_complex(BMO.encode_complex(1.5 - 2im)) == 1.5 - 2im
    end

    @testset "Assets" begin
        ctx = BMO.StorageContext(v"0.0.1")
        bytes = rand(UInt8, 100)
        p1 = BMO.write_asset!(ctx, "mesh", "bin", bytes)
        p2 = BMO.write_asset!(ctx, "other", "bin", copy(bytes))
        @test p1 == p2
        @test startswith(p1, "assets/mesh-")
        @test length(ctx.assets) == 1
        @test_throws ArgumentError BMO.write_asset!(ctx, "a/b", "bin", rand(UInt8, 3))
        @test_throws r"asset \"assets/missing.bin\" is missing" BMO.read_asset(BMO.StorageContext(v"0.0.1"), "assets/missing.bin")

        # Identical large arrays of two objects are stored once
        data = rand(3, 1000)
        path = tmpfile()
        BMO.save_setup(path; systems = [System([StoreArrayObj(data), StoreArrayObj(copy(data))])])
        setup = BMO.load_setup(path)
        loaded = setup.systems[1].objects
        @test loaded[1].data == data && loaded[2].data == data
        reader = BMO.ZipReader(read(path))
        @test count(n -> startswith(n, "assets/"), BMO.zip_names(reader)) == 1
    end

    @testset "Systems, pairs and sharing" begin
        shared = StoreObj(BMO.Point3(0.1 + 0.2, 1 / 3, 1e-300))
        a, b = StoreObj(BMO.Point3(1.0, 2, 3)), StoreObj(BMO.Point3(-1.0, 0, 0))
        sys1, sys2 = System([a, shared]), System([shared, b])
        src1, src2 = StoreSrc(1e-6), StoreSrc(532e-9)
        extra_sys, extra_src = System([StoreObj(BMO.Point3(0.0, 0, 0))]), StoreSrc(1.55e-6)
        names = IdDict{Any, String}(sys1 => "bench", shared => "shared", src1 => "laser")
        metadata = Dict("gui" => Dict("camera" => Dict("eye" => [1.0, 2.0, 3.0]), "labels" => Dict("shared" => "M1")))
        path = tmpfile()
        BMO.save_setup(path, sys1 => src1, sys2 => src1, sys2 => src2;
            systems = [extra_sys], sources = [extra_src], names, metadata)

        setup = BMO.load_setup(path)
        @test setup.version == pkgversion(BeamletOptics)
        @test length(setup.pairs) == 3
        @test length(setup.systems) == 3
        @test length(setup.sources) == 3
        (s1, l1), (s2, l2), (s2b, l3) = setup.pairs
        @test s2 === s2b
        @test l1 === l2
        @test l1 !== l3
        @test l1.λ == 1e-6 && l3.λ == 532e-9
        # The shared object is the same instance in both systems
        @test s1.objects[2] === s2.objects[1]
        @test setup.names["shared"] === s1.objects[2]
        @test setup.names["bench"] === s1
        @test setup.names["laser"] === l1
        @test s1.objects[2].pos == shared.pos
        @test s1.objects[1].pos == a.pos && s2.objects[2].pos == b.pos
        @test setup.metadata == metadata
        @test setup.sources[3].λ == 1.55e-6
        @test setup.systems[3].objects[1].pos == BMO.Point3(0.0, 0, 0)

        @test_throws r"duplicate name" BMO.save_setup(tmpfile(), sys1 => src1; names = IdDict(a => "x", b => "x"))
        @test_throws r"cannot be stored" BMO.save_setup(tmpfile(); systems = [System([Unregistered()])])
    end

    @testset "Several sources per system" begin
        sys1, sys2 = System([StoreObj(BMO.Point3(1.0, 0, 0))]), System([StoreObj(BMO.Point3(2.0, 0, 0))])
        src1, src2, src3 = StoreSrc(1e-6), StoreSrc(532e-9), StoreSrc(1.55e-6)
        path = tmpfile()
        BMO.save_setup(path, sys1 => (src1, src2), sys2 => [src2, src3], sys2 => (src1,), sys1 => src3)
        setup = BMO.load_setup(path)
        (s1, t1), (s2, t2), (s2b, l3), (s1b, l4) = setup.pairs
        @test s1 === s1b && s2 === s2b
        @test t1 isa Tuple && length(t1) == 2 && t2 isa Tuple && length(t2) == 2
        @test [s.λ for s in t1] == [1e-6, 532e-9]
        @test [s.λ for s in t2] == [532e-9, 1.55e-6]
        # Sources are shared across pairs, single sources are returned unwrapped
        @test t1[2] === t2[1]
        @test l3 === t1[1]
        @test l4 === t2[2]
        @test length(setup.sources) == 3

        @test_throws r"at least one source" BMO.save_setup(tmpfile(), sys1 => ())
    end

    @testset "Object groups" begin
        a, b, c = StoreObj(BMO.Point3(1.0, 0, 0)), StoreObj(BMO.Point3(0.0, 1, 0)), StoreObj(BMO.Point3(0.0, 0, 1))
        inner = ObjectGroup([b, c])
        inner.center = BMO.Point3(0.5, 0.5, 0.5)
        inner.dir = BMO.SMatrix{3, 3}(BMO.rotate3d([0, 0, 1], 0.3))
        outer = ObjectGroup([a, inner])
        outer.center = BMO.Point3(-1.0, 2.0, 0.25)
        outer.dir = BMO.SMatrix{3, 3}(BMO.rotate3d([1, 1, 0] / sqrt(2), 1.1))
        path = tmpfile()
        BMO.save_setup(path; systems = [System([outer])], names = IdDict{Any, String}(c => "c"))
        setup = BMO.load_setup(path)
        louter = setup.systems[1].objects[1]
        @test louter isa ObjectGroup
        @test louter.center == outer.center && louter.dir == outer.dir
        linner = louter.objects[2]
        @test linner.center == inner.center && linner.dir == inner.dir
        @test louter.objects[1].pos == a.pos
        @test [o.pos for o in linner.objects] == [b.pos, c.pos]
        @test setup.names["c"] === linner.objects[2]

        # StaticSystems are stored flattened
        path = tmpfile()
        BMO.save_setup(path; systems = [StaticSystem([a, b])])
        @test [o.pos for o in BMO.load_setup(path).systems[1].objects] == [a.pos, b.pos]
    end

    @testset "Load errors" begin
        path = tmpfile()
        BMO.save_setup(path, System([StoreObj(BMO.Point3(0.0, 0, 0))]) => StoreSrc(1e-6))
        v = pkgversion(BeamletOptics)
        other = v.major == 0 ? VersionNumber(0, v.minor + 1, 0) : VersionNumber(v.major + 1, 0, 0)
        @test_throws r"use|release" BMO.load_setup(path; _loader = other)
        @test !isnothing(BMO.load_setup(path; _loader = v))

        delete!(BMO.STORAGE_TYPES, "StoreSrc")
        try
            @test_throws r"unknown type \"StoreSrc\" in setup written by BeamletOptics" BMO.load_setup(path)
        finally
            BMO.register_storage_type!(StoreSrc, "StoreSrc")
        end
        @test_throws r"already registered" BMO.register_storage_type!(StoreObj, "StoreSrc")

        notzip = joinpath(mktempdir(), "x.bmo")
        write(notzip, "hello")
        @test_throws r"not a BeamletOptics setup file" BMO.load_setup(notzip)
        @test_throws r"does not exist" BMO.load_setup(joinpath(mktempdir(), "none.bmo"))
    end
end

end # module
