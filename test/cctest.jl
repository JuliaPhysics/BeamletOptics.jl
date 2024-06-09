using GLMakie, SCDI, GeometryBasics

## define components and positions
bs = SCDI.ThinBeamSplitter(40e-3);
pd_1 = SCDI.Photodetector(10e-3, 250);
pd_2 = SCDI.Photodetector(10e-3, 250);

SCDI.zrotate3d!(bs, deg2rad(45))
SCDI.translate3d!(pd_1, [0, 0.2, 0])
SCDI.zrotate3d!(pd_1, deg2rad(180))

SCDI.translate3d!(pd_2, [0.2, 0, 0])
SCDI.zrotate3d!(pd_2, deg2rad(90))

# define system and beams -> solve
system = SCDI.StaticSystem([bs, pd_1, pd_2]);

#SCDI.zrotate3d!(bs, deg2rad(0.025)) # ENABLE to create bug

phis = LinRange(0, 2pi, 50)
p1 = similar(phis)
p2 = similar(phis)

for (i, phi) in enumerate(phis)
    global l1 = SCDI.GaussianBeamlet(SCDI.Ray([0, -0.1, 0], [0, 1., 0]), 1064e-9, .5e-3, P0=0.5);
    global l2 = SCDI.GaussianBeamlet(SCDI.Ray([-0.1, 0, 0], [1., 0, 0]), 1064e-9, .5e-3, P0=0.5);
    l1.E0 *= exp(im*phi)
    SCDI.reset_photodetector!(pd_1)
    SCDI.reset_photodetector!(pd_2)
    SCDI.solve_system!(system, l1)
    SCDI.solve_system!(system, l2)
    p1[i] = SCDI.optical_power(pd_1)
    p2[i] = SCDI.optical_power(pd_2)
end

##
pfig, pax = lines(phis, p1)
lines!(pax, phis, p2)
lines!(pax, phis, p1.+p2)
pax2 = Axis(pfig[2,1])
lines!(pax2, phis, p1 .+ p2 .- 1)

pfig

## plot
fig = Figure()
ax = LScene(fig[1,1:4])
SCDI.render_system!(ax, system)
SCDI.render_beam!(ax, l1)
SCDI.render_beam!(ax, l2)

SCDI.render_object_normals!(ax, pd_1.shape)
SCDI.render_object_normals!(ax, pd_2.shape)
SCDI.render_object_normals!(ax, bs.shape)

heat1 = Axis3(fig[2, 1], xlabel="x [m]", ylabel="y [m]", title="P = $(SCDI.optical_power(pd_1)) W")
hm1 = surface!(heat1, pd_1.x, pd_1.y, SCDI.intensity(pd_1))
cb1 = Colorbar(fig[2, 2], hm1, label="Intensity [W/m²]")

heat = Axis3(fig[2, 3], xlabel="x [m]", ylabel="y [m]", title="P = $(SCDI.optical_power(pd_2)) W")
hm = surface!(heat, pd_2.x, pd_2.y, SCDI.intensity(pd_2))
cb = Colorbar(fig[2, 4], hm, label="Intensity [W/m²]")

@info 1 - (SCDI.optical_power(pd_1) + SCDI.optical_power(pd_2))

fig


##

function foo(root)
    for node in StatelessBFS(root)
        @info "Hi"
    end

    for node in Leaves(root)
        @info "Yo"
    end
end

##
using AbstractTrees

struct MyNode
    value::Int
    children::Vector{MyNode}
end

# Implement the AbstractTrees interface
AbstractTrees.children(node::MyNode) = node.children

# Create a sample tree
leaf1 = MyNode(1, [])
leaf2 = MyNode(2, [])
root = MyNode(0, [leaf1, leaf2])

function stateful_bfs(root::MyNode)
    queue = MyNode[root]  # Initialize the queue with the root node
    while !isempty(queue)
        node = popfirst!(queue)  # Dequeue the first node
        println(node.value)  # Process the node (print its value in this case)
        append!(queue, children(node))  # Enqueue the children
    end
end

stateful_bfs(root)
