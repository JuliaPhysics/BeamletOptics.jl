using GLMakie, BeamletOptics
GLMakie.activate!(; ssao=true) #hide
const BMO = BeamletOptics
const mm = 1e-3

λ = 632.8e-9             # HeNe wavelength
Δx = 100mm               # lateral offset between laser and optical axis

function mounted_mirror()
    mirror = RoundPlanoMirror(BMO.inch, 6mm)
    return ObjectGroup([BMO.KM100CPMount(), mirror])
end

m1 = mounted_mirror()
m2 = mounted_mirror()
zrotate3d!(m1, deg2rad(45));   translate3d!(m1, [0, 100mm, 0])
zrotate3d!(m2, deg2rad(-135)); translate3d!(m2, [Δx, 100mm, 0])

function optical_table!(ax; xs=(-75mm, 175mm), ys=(-50mm, 500mm), z_top=-81.8mm, pitch=25mm)
    mesh!(ax, Rect3f(Vec3f(xs[1], ys[1], z_top - 12mm), Vec3f(xs[2] - xs[1], ys[2] - ys[1], 12mm)), color=:grey85)
    holes = [Point3f(x, y, z_top + 0.2mm) for x in xs[1]+pitch/2:pitch:xs[2], y in ys[1]+pitch/2:pitch:ys[2]]
    meshscatter!(ax, vec(holes), markersize=Vec3f(3mm, 3mm, 0.1mm), color=:grey45)   # flattened spheres: Ø6 mm holes
end

system = System([m1, m2])
beam = Beam(Ray([0, 0, 0], [0, 1.0, 0], λ))
solve_system!(system, beam)

fig = Figure(size=(600, 400))
ax = LScene(fig[1,1], show_axis=false)
optical_table!(ax)
render!(ax, system)
render!(ax, beam, color=:red, flen=0.3, show_pos=true)
render_lcs!(ax; scale=8, show_labels=true)
set_view(ax, [212mm, -106mm, 138mm], [42mm, 87mm, -46mm], [0, 0, 1]) #hide
fig

using Base64
# The mounting error of M1 is encoded on purpose. Try the exercise below before decoding it!
θx, θz = deg2rad.(parse.(Float64, split(String(base64decode("LTAuMjMgMC4zNw==")))))
xrotate3d!(m1, θx); zrotate3d!(m1, θz)

card = Detector(BMO.inch, false)       # alignment card: records hits, lets the beam pass
translate3d!(card, [Δx, 300mm, 0])
system = System([m1, m2, card])

function check_alignment()
    empty!(card)
    solve_system!(system, Beam(Ray([0, 0, 0], [0, 1.0, 0], λ)))
    offset = spot_diagram(card)[1]     # [x, z] on the card
    println("offset on card: x = ", round(offset[1] / mm, digits=4) + 0, " mm, z = ", round(offset[2] / mm, digits=4) + 0, " mm")
    return offset
end

offset = check_alignment()

function plot_alignment()
    beam = Beam(Ray([0, 0, 0], [0, 1.0, 0], λ))
    solve_system!(system, beam)
    fig = Figure(size=(600, 400))
    ax = LScene(fig[1,1], show_axis=false)
    optical_table!(ax)
    render!(ax, system)
    render!(ax, beam, color=:red, flen=0.3, show_pos=true)
    spot_ax = Axis(fig[1,2], aspect=1, xlabel="x [mm]", ylabel="z [mm]",
        limits=(-12.7, 12.7, -12.7, 12.7), title="Card")
    scatter!(spot_ax, [0.0], [0.0], color=:black, marker=:xcross, markersize=16)
    pts = spot_diagram(card)
    scatter!(spot_ax, [p[1]/mm for p in pts], [p[2]/mm for p in pts], color=:red)
    set_view(ax, [157mm, 465mm, 106mm], [77mm, 208mm, -64mm], [0, 0, 1]) #hide
    return fig
end

fig2 = plot_alignment()
fig2

L = 300mm                              # M1 → M2 → card
# Your turn, e.g.:
# zrotate3d!(m1, ...)
# xrotate3d!(m1, ...)
# check_alignment()


# SPOILER: solution below. Try the exercise first!
zrotate3d!(m1, deg2rad(-0.37))
xrotate3d!(m1, deg2rad(0.23))
check_alignment()

fig3 = plot_alignment()
fig3
