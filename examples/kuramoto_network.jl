#=
## Integration with `NetworkDynamics.jl`
In this example we model a Network based on the [Kuroamoto model](https://en.wikipedia.org/wiki/Kuramoto_model) using the [`NetworkDynamics.jl`](https://github.com/FHell/NetworkDynamics.jl) package.
=#

using NetworkDynamics
using LightGraphs
using IOSystems
using ModelingToolkit
using DifferentialEquations
using Plots

#=
In the Kuramoto model each vertex $i$ has it's own intrinsic angular frequency $\omega_i$.
The angle of frequency is given by
```math
\dot\phi = \omega_i + \sum e
```
where $\sum e$ sums over all the connected edges.
=#
@parameters t ω edgesum(t)
@variables ϕ(t)
@derivatives D'~t

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

## allocation free oop aggregator. might be more difficult for more-dimensional systems
## unfortunately there are no vector variables in MDK and we can't model the aggregators
## as IOSystems.
function sumedges(edges)
    r = 0.0
    for e in edges
        r += e[1]
    end
    return r
end

vertex = generate_ode_vertex(kmvert, sumedges)
nothing #hide

#=
The edges are defined as
```math
e_{i,j} = K \cdot \mathrm{sin}(\phi_i - \phi_j)
```
For the (static) edges we can't used IOSystem function building in the moment.
Well, we could but than we'd get an `ODEEdge` with zero-massmatrix.
=#
@parameters v₁(t) v₂(t) K
edge_ip = build_function([K*sin(v₁[1] - v₂[1])], v₁, v₂, K, t, expression=Val{false})[2]
edge = StaticEdge(f! = edge_ip, dim = 1, coupling=:antisymmetric)

nothing #hide
#=
Lets runt the [same example](https://fhell.github.io/NetworkDynamics.jl/dev/heterogeneous_system/#Heterogenous-parameters) as in `NetworkDynamics.jl`
=#
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
sol = solve(prob)
plot(sol, ylabel="ϕ")
