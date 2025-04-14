function render!(
        ax::_RenderEnv,
        mesh::BMO.AbstractMesh;
        show_normals::Bool=false,
        show_normals_length::Real=0.01,
        kwargs...
    )
    # Plot mesh
    mesh!(ax, BMO.vertices(mesh), BMO.faces(mesh); kwargs...)
    # Plot mesh normals
    if show_normals
        for fID in 1:size(BMO.faces(mesh))[1]
            nml = BMO.normal3d(mesh, fID)
            pos = @view BMO.vertices(mesh)[BMO.faces(mesh)[fID, 1], :]
            vec = pos + show_normals_length * nml
            lines!(
                ax,
                [pos[1], vec[1]],
                [pos[2], vec[2]],
                [pos[3], vec[3]];
                color=:blue,
            )
        end
    end
    return nothing
end