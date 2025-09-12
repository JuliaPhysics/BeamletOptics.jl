using CairoMakie, BeamletOptics

CairoMakie.activate!()

function generate_fresnel_plot(n1, n2)
    n = n2/n1
    θ = deg2rad.(0:.01:90)
    rs, rp, ts, tp = BeamletOptics.fresnel_coefficients(θ, n)
    
    coeff = Figure(size=(600,300))
    ax1 = Axis(coeff[1,1], xlabel="Angle of incidence θ [°]", ylabel="Real.", title="n₁ = $n1, n₂ = $n2")
    ax2 = Axis(coeff[1,1], xlabel="Angle of incidence θ [°]", ylabel="Imag.", yaxisposition = :right)
    
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

    return coeff
end

##
n1 = 1.0
n2 = 1.5
fig = generate_fresnel_plot(n1, n2)
save("vac_to_glass.png", fig, px_per_unit=4)

n1 = 1.5
n2 = 1.0
fig = generate_fresnel_plot(n1, n2)
save("glass_to_vac.png", fig, px_per_unit=4)
