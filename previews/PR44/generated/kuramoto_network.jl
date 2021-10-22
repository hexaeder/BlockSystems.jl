using NetworkDynamics
using LightGraphs
using BlockSystems
using ModelingToolkit
using OrdinaryDiffEq
using Plots

@parameters t ω edgesum(t)
@variables ϕ(t)
D = Differential(t)

kmvert = IOBlock([D(ϕ) ~ ω + edgesum],
                 [edgesum],
                 [ϕ],
                 name=:kmvert)

function generate_ode_vertex(iob::IOBlock, aggregator; kwargs...)
    gen = generate_io_function(iob; kwargs...)
    f = (du, u, edges, p, t) -> gen.f_ip(du, u, aggregator(edges), p, t)
    vars = ModelingToolkit.tosymbol.(gen.states)
    ODEVertex(f! = f, dim = length(vars), sym = vars)
end

# allocation free oop aggregator. might be more difficult for more-dimensional systems
# unfortunately there are no vector variables in MDK and we can't model the aggregators
# as an IOSystem.
function sumedges(edges)
    r = 0.0
    for e in edges
        r += e[1]
    end
    return r
end

vertex = generate_ode_vertex(kmvert, sumedges)
nothing #hide

@parameters v₁(t) v₂(t) K
edge_ip = build_function([K*sin(v₁[1] - v₂[1])], [v₁], [v₂], K, t, expression=Val{false})[2]
edge = StaticEdge(f! = edge_ip, dim = 1, coupling=:antisymmetric)

nothing #hide

N = 8
g = watts_strogatz(N,2,0) # ring network
nd = network_dynamics(vertex, edge, g)

ω = collect(1:N)./N
ω .-= sum(ω)/N
K = 3.
p = (ω, K); # p[1] vertex parameters, p[2] edge parameters

x0 = collect(1:N)./N
x0 .-= sum(x0)./N


tspan = (0., 10.)
prob = ODEProblem(nd, x0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, ylabel="ϕ")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

