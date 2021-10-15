using LightGraphs
using BlockSystems
using ModelingToolkit
using OrdinaryDiffEq
using Plots

function gen_edge_block(name)
    @parameters t src(t) dst(t) K
    @variables out(t)
    IOBlock([out ~ K*sin(src-dst)], [src, dst], [out], name=Symbol(name))
end

# append subscript, `sybscript(:foo, 2) => :foo₂`
subscript(s, i) = Symbol(s, Char(0x02080 + i))

function gen_vertex_block(n_edges, name)
    @parameters t ω
    @variables ϕ(t)
    D = Differential(t)
    # the way array variables work changed. This is a hack to retrieve the old behvaviour of
    # this closely mimics the old @parameters edge[1:n_edges](t)
    edge = Num[]
    for i in 1:n_edges
        symname = subscript(:edge, i)
        append!(edge, @parameters $symname(t))
    end

    IOBlock([D(ϕ) ~ ω + (+)(edge...)],
            [edge...],
            [ϕ],
            name=Symbol(name))
end
nothing #hide

N = 8
g = SimpleDiGraph(watts_strogatz(N,2,0)) # ring network
nothing #hide

edgelist = [(i=i, src=e.src, dst=e.dst, block=gen_edge_block("e_$(e.src)_$(e.dst)"))
            for (i, e) in enumerate(edges(g))]
edge_blocks = [e.block for e in edgelist]
nothing #hide

vert_blocks = IOBlock[]
connections = Pair[]

for i in vertices(g)
    # collect the incoming edges for each node
    edges = filter(e -> e.dst == i, edgelist)

    node = gen_vertex_block(length(edges), "v$i")
    push!(vert_blocks, node)

    # each node has the open inputs edge₁, edge₂, ...
    # we need to connect the outputs of the edge-blocks to the
    # inputs of the node like edge_j_to_1.out => node.edge₁
    for (i, edge) in enumerate(edges)
        node_input_i = getproperty(node, subscript(:edge, i))
        push!(connections, edge.block.out => node_input_i)
    end
end

for edge in edgelist
    push!(connections, vert_blocks[edge.src].ϕ => edge.block.src)
    push!(connections, vert_blocks[edge.dst].ϕ => edge.block.dst)
end

outputs = [block.ϕ for block in vert_blocks]

network = IOSystem(connections, vcat(vert_blocks, edge_blocks), outputs=outputs)

networkblock = connect_system(network, verbose=false)
nothing #hide

gen = generate_io_function(networkblock,
                           f_states=[v.ϕ for v in vert_blocks],
                           f_params=vcat([v.ω for v in vert_blocks],
                                         [e.K for e in edge_blocks]),
                           warn=false);
nothing #hide

odefun(du, u, p, t) = gen.f_ip(du, u, (), p, t)
nothing #hide

ω = collect(1:N)./N
ω .-= sum(ω)/N
K = [3.0 for i in edge_blocks]
p = (ω..., K...)

x0 = collect(1:N)./N
x0 .-= sum(x0)./N
nothing #hide

tspan = (0., 10.)
prob = ODEProblem(odefun, x0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, ylabel="ϕ")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

