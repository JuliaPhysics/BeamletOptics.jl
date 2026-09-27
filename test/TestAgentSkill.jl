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

    @testset "install_agent_skill copies the shipped skill" begin
        mktempdir() do dest
            skill_md = joinpath(dest, "beamletoptics", "SKILL.md")
            original = read(joinpath(SKILL_DIR, "SKILL.md"), String)

            path = @test_logs (:info, r"Installed BeamletOptics agent skill") BeamletOptics.install_agent_skill(dest)
            @test path == joinpath(dest, "beamletoptics")
            @test read(skill_md, String) == original
            @test isfile(joinpath(path, "scripts", "api_lookup.jl"))
            @test isfile(joinpath(path, "templates", "01_singlet_spot_diagram.jl"))

            # a reinstall replaces an edited, read-only copy (as left behind by a read-only Pkg install)
            write(skill_md, "edited")
            for (root, _, files) in walkdir(path)
                foreach(f -> chmod(joinpath(root, f), 0o444), files)
            end
            @test_logs (:info,) BeamletOptics.install_agent_skill(dest)
            @test read(skill_md, String) == original
            @test uperm(skill_md) & 0x02 != 0   # writable again
        end
    end
end

end
