# Writes the fixture `michelson.bmo` that TestStorageE2E.jl loads.
# Regenerate it after a breaking release (new major version, or new minor version for 0.x):
#   julia --project=test test/Storage/fixtures/make_fixture.jl
using BeamletOptics

const BMO = BeamletOptics

l_0 = 0.1
m1 = SquarePlanoMirror2D(BMO.inch)
m2 = SquarePlanoMirror2D(BMO.inch)
bs = ThinBeamsplitter(BMO.inch, reflectance = 0.5)
pd = Detector(BMO.inch / 5)
translate3d!(m1, [l_0, 0, 0])
translate3d!(m2, [0, l_0, 0])
translate3d!(pd, [-l_0, 0, 0])
zrotate3d!(bs, deg2rad(45))
zrotate3d!(m1, deg2rad(90))
zrotate3d!(pd, deg2rad(90))
beam = GaussianBeamlet([0, -l_0, 0], [0, 1.0, 0], 635e-9, 1e-4, P0 = 5e-3)

names = IdDict{Any, String}(m1 => "m1", m2 => "m2", bs => "bs", pd => "pd", beam => "laser")
save_setup(joinpath(@__DIR__, "michelson.bmo"), System([m1, m2, bs, pd]) => beam; names,
    metadata = Dict("description" => "Michelson interferometer fixture"))
