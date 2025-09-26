# order or inclusion matters!
include(joinpath(@__DIR__, "TestUtils.jl"))
include(joinpath(@__DIR__, "TestAbstractTypes.jl"))

# Test ray and beam types
include(joinpath(@__DIR__, "TestRays.jl"))
include(joinpath(@__DIR__, "TestPolarizedRays.jl"))
include(joinpath(@__DIR__, "TestBeams.jl"))
include(joinpath(@__DIR__, "TestBeamGroups.jl"))
include(joinpath(@__DIR__, "TestGaussianBeamlet.jl"))

# Test geometry representation
include(joinpath(@__DIR__, "Geometry", "TestMesh.jl"))
include(joinpath(@__DIR__, "Geometry", "TestSDFs.jl"))

# Test system and object containers
include(joinpath(@__DIR__, "TestSystem.jl"))
include(joinpath(@__DIR__, "TestObjectGroups.jl"))

# Test lens models
include(joinpath(@__DIR__, "Lenses", "TestSphericalLenses.jl"))
include(joinpath(@__DIR__, "Lenses", "TestSurfaces.jl"))
include(joinpath(@__DIR__, "Lenses", "TestAsphericalLenses.jl"))
include(joinpath(@__DIR__, "Lenses", "TestCylindricalLenses.jl"))

# Test component models
include(joinpath(@__DIR__, "Components", "TestDummies.jl"))
include(joinpath(@__DIR__, "Components", "TestDetector.jl"))
include(joinpath(@__DIR__, "Components", "TestBeamsplitters.jl"))
include(joinpath(@__DIR__, "Components", "TestPolarizers.jl"))

# Test end-to-end models
include(joinpath(@__DIR__, "E2E", "TestDoubleGaussLens.jl"))
include(joinpath(@__DIR__, "E2E", "TestMichelson.jl"))
include(joinpath(@__DIR__, "E2E", "TestMachZehnder.jl"))

# Test regressions
include(joinpath(@__DIR__, "TestBugFixes.jl"))

include(joinpath(@__DIR__, "TestMisc.jl"))
