using LightGraphs
using BlockSystems
using ModelingToolkit
using OrdinaryDiffEq
using Plots

function gen_edge_block(name)
    @parameters t src(t) dst(t) K
    @variables o(t)
    IOBlock([ o ~ K*sin(src-dst) ], [src,dst], [o], name=Symbol(name))
end

function gen_vertex_block(n_in, n_out, name)
    @parameters t ω edgesum(t)
    @parameters in_edge[1:n_in](t)
    @parameters out_edge[1:n_out](t)
    @variables ϕ(t)
    D = Differential(t)

    rhs = ω
    if n_in > 0
        rhs += (+)(in_edge...)
    end
    if n_out > 0
        rhs -= (+)(out_edge...)
    end

    IOBlock([D(ϕ) ~ rhs],
            [in_edge..., out_edge...],
            [ϕ],
            name=Symbol(name))
end
nothing #hide

N = 8
g = watts_strogatz(N,2,0) # ring network
nothing #hide

edgelist = [(i=i, src=e.src, dst=e.dst, block=gen_edge_block("$(e.src)to$(e.dst)")) for (i, e) in enumerate(edges(g))]
edge_blocks = [e.block for e in edgelist]
nothing #hide

vert_blocks = IOBlock[]
connections = Pair[]

for i in vertices(g)
    out_edges = filter(e->e.src == i, edgelist)
    in_edges = filter(e->e.dst == i, edgelist)
    block = gen_vertex_block(length(in_edges), length(out_edges), "v$i")
    push!(vert_blocks, block)

    for (i, edge) in enumerate(out_edges)
        input = getproperty(block, Symbol("out_edge", Char(0x02080 + i)))
        push!(connections, edge.block.o => input)
    end

    for (i, edge) in enumerate(in_edges)
        input = getproperty(block, Symbol("in_edge", Char(0x02080 + i)))
        push!(connections, edge.block.o => input)
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
                                         [e.K for e in edge_blocks]))
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

