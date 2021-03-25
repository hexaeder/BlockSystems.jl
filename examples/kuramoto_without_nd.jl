#=
## Network Dynamics without NetworkDynamics.jl
=#
## In this example we want to explore the same problem as in [Integration with `NetworkDynamics.jl`](@ref).
## But this time without `NetworkDynamics.jl` ...

using LightGraphs
using BlockSystems
using ModelingToolkit
using OrdinaryDiffEq
using Plots

#=
The goal is to generate IOBlocks for edges and vertices based on a given graph.
In order to do so we have to make the vertex blocks a bit special: since MDK
does not support vector inputs we need a special IOBlock which depeends on the number
of connected edges.

In our model we look at a directed graph. Each edge is represented by a function on such edge.
Each edge funciton sees the values of the connected nodes.

The nodes are also modeld as functions. Each node sees all the *incoming* edge values.

```
     e₁₂ = f(1,2)
     .---→---.
   (1)       (2)
     ˙---←---˙
     e₂₁ = f(2,1)
```
=#
function gen_edge_block(name)
    @parameters t src(t) dst(t) K
    @variables o(t)
    IOBlock([o ~ K*sin(src-dst)], [src, dst], [o], name=Symbol(name))
end

function gen_vertex_block(n_edges, name)
    @parameters t ω
    @parameters edge[1:n_edges](t)
    @variables ϕ(t)
    D = Differential(t)

    IOBlock([D(ϕ) ~ ω + (+)(edge...)],
            [edge...],
            [ϕ],
            name=Symbol(name))
end
nothing #hide

#=
As a graph we chose a ring network of

=#
N = 8
g = SimpleDiGraph(watts_strogatz(N,2,0)) # ring network
nothing #hide

# First we generate a list of all edge-blocks because the don't depend on the vertices.
edgelist = [(i=i, src=e.src, dst=e.dst, block=gen_edge_block("e_$(e.src)_$(e.dst)")) for (i, e) in enumerate(edges(g))]
edge_blocks = [e.block for e in edgelist]
nothing #hide

#=
With the edges we can generate vert blocks based on their number of out and in edges.
We can also create the connections
```
e1to2.src => v1.out_edge₁
e1to2.dst => v2.in_edge₁
```
and so forth.
=#
vert_blocks = IOBlock[]
connections = Pair[]

for i in vertices(g)
    # collect the incoming edges for each node
    edges = filter(e -> e.dst == i, edgelist)

    node = gen_vertex_block(length(edges), "v$i")
    push!(vert_blocks, node)

    # each node has the open inputs edge₁, edge₂, ...
    # we need to connect the ouputs of those edge-blocks to the
    # inputs of the node like edge1.o => node.edge₁
    for (i, edge) in enumerate(edges)
        node_input_i = getproperty(node, Symbol("edge", Char(0x02080 + i)))
        push!(connections, edge.block.o => node_input_i)
    end
end

# Once the vertices are generated we can plug the edges src and dst to the output of the corresponding vertex.
for edge in edgelist
    push!(connections, vert_blocks[edge.src].ϕ => edge.block.src)
    push!(connections, vert_blocks[edge.dst].ϕ => edge.block.dst)
end

#=
We want the `connect_system` to get rid of the algebraic states for the edges. Therefore
we have to provide a list of `outputs` which only contains the outputs of the vertices. By doing
so the edge outputs will become `istates` of the `IOSystem` and upon connection might be reduced.
=#
outputs = [block.ϕ for block in vert_blocks]

network = IOSystem(connections, vcat(vert_blocks, edge_blocks), outputs=outputs)

networkblock = connect_system(network, verbose=false)
nothing #hide

# As the output shows the system has be reduced to just N equations.
# Well now we can generate the functions...
gen = generate_io_function(networkblock,
                           f_states=[v.ϕ for v in vert_blocks],
                           f_params=vcat([v.ω for v in vert_blocks],
                                         [e.K for e in edge_blocks]))
nothing #hide

# ... enclose the `f_ip` to get rid of the empty `inputs` field...
odefun(du, u, p, t) = gen.f_ip(du, u, (), p, t)
nothing #hide

# ... set the starting conditions ...
ω = collect(1:N)./N
ω .-= sum(ω)/N
K = [3.0 for i in edge_blocks]
p = (ω..., K...)

x0 = collect(1:N)./N
x0 .-= sum(x0)./N
nothing #hide

#src using BenchmarkTools
#src @btime $odefun($x0, $x0, $p, 0.0)

# ... and solve the system!
tspan = (0., 10.)
prob = ODEProblem(odefun, x0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, ylabel="ϕ")
