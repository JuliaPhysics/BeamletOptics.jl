function render!(ax::_RenderEnv, mesh::BMO.AbstractMesh; kwargs...)
    mesh!(ax, BMO.vertices(mesh), BMO.faces(mesh); kwargs...)
    return nothing
end