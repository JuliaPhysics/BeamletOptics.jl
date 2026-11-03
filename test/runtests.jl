# Every file in test/ is a self-contained test module, run in parallel worker processes.
# Usage: julia --project=test test/runtests.jl [--jobs=N] [--list] [TESTS...]
using BeamletOptics
using ParallelTestRunner

# MUST run in a fresh worker: tests the errors of render! before the BMO Makie ext is loaded
const makie_free_tests = ("Rendering/TestRenderErrors",)

worker_of(test::String) = test in makie_free_tests ? addworker() : nothing

runtests(BeamletOptics, ARGS; testsuite = find_tests(@__DIR__), test_worker = worker_of)
