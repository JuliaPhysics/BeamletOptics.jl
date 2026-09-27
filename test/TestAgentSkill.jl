module TestAgentSkill

using BeamletOptics
using Test

const SKILL_DIR = joinpath(@__DIR__, "..", "skills", "beamletoptics")

@testset "Agent skill" begin
    @testset "API.md export table matches names(BeamletOptics)" begin
        api = read(joinpath(SKILL_DIR, "API.md"), String)
        table = split(split(api, "## Exported names")[2], "\n## ")[1]
        rows = filter(l -> startswith(l, "|") && !occursin(r"^\|\s*(Category|-)", l),
            split(table, '\n'))
        documented = Set{Symbol}()
        for row in rows, m in eachmatch(r"`([A-Za-z_][A-Za-z0-9_!]*)`", split(row, '|')[3])
            push!(documented, Symbol(m.captures[1]))
        end
        exported = Set(filter(!=(:BeamletOptics), names(BeamletOptics)))
        @test isempty(setdiff(documented, exported))  # documented but not exported
        @test isempty(setdiff(exported, documented))  # exported but not documented
    end

    @testset "SKILL.md targets the current minor version" begin
        skill = read(joinpath(SKILL_DIR, "SKILL.md"), String)
        m = match(r"beamletoptics-version:\s*\"(\d+)\.(\d+)\"", skill)
        @test m !== nothing
        v = pkgversion(BeamletOptics)
        @test parse.(Int, m.captures) == [v.major, v.minor]
    end
end

end
