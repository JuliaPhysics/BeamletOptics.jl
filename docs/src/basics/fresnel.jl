function plot_and_save_fresnel_coeffs(n1, n2; save_fig::Bool=false, path::String=@__DIR__)
    n = n2/n1
    θ = deg2rad.(0:.01:90)
    rs, rp, ts, tp = SCDI.fresnel_coefficients(θ, n)

    coeff = Figure(size=(500,300))
    ax1 = Axis(coeff[1,1], xlabel="θ [°]", ylabel="Real.", title="n₁ = $n1, n₂ = $n2")
    ax2 = Axis(coeff[1,1], xlabel="θ [°]", ylabel="Imag.", yaxisposition = :right)

    hidespines!(ax2)
    hidexdecorations!(ax2)
    # linkxaxes!(ax1, ax2)
    linkyaxes!(ax1, ax2)

    xlims!(ax1, minimum(rad2deg.(θ)), maximum(rad2deg.(θ)))
    xlims!(ax2, minimum(rad2deg.(θ)), maximum(rad2deg.(θ)))

    lines!(ax1, rad2deg.(θ), real.(rs), color=:orange, label="rs")
    lines!(ax2, rad2deg.(θ), imag.(rs), color=:orange, linestyle=:dashdot)

    lines!(ax1, rad2deg.(θ), real.(rp), color=:red, label="rp")
    lines!(ax2, rad2deg.(θ), imag.(rp), color=:red, linestyle=:dashdot)

    lines!(ax1, rad2deg.(θ), real.(ts), color=:dodgerblue, label="ts")
    lines!(ax2, rad2deg.(θ), imag.(ts), color=:dodgerblue, linestyle=:dashdot)

    lines!(ax1, rad2deg.(θ), real.(tp), color=:blue, label="tp")
    lines!(ax2, rad2deg.(θ), imag.(tp), color=:blue, linestyle=:dashdot)

    # axislegend(ax1, ax1, position=:lt)
    coeff[1, 2] = Legend(coeff, ax1, "Fresnel (real)", framevisible = false)

    resize_to_layout!(coeff)
    
    if save_fig
        save("fresnel_coeffs_$(n1)_$(n2).png", coeff, px_per_unit=4)
    end
    return coeff
end