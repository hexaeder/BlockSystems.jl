#=
## Generate `PowerDynamics.jl`-Node with `BlockSystems`
We start by defining some general stuff...
=#
using BlockSystems
using ModelingToolkit

@parameters t
D = Differential(t)
nothing #hide

#=
### Defining the common Blocks
We need some common parts such as filters or integrators. The actual Blocks will
be defined by using these 'blueprints'. Stuff like this should go to a library of
components at some point.

=#
# #### low pass filter
@parameters τ input(t)
@variables filtered(t)

lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)],
              [input], [filtered])

# #### integrator
@parameters x(t)
@variables int(t)

integrator = IOBlock([D(int) ~ x], [x], [int])

# #### voltage source

@parameters ω(t) v(t) τ
@variables u_i(t) u_r(t)

voltage_source = IOBlock([D(u_r) ~ -ω * u_i + (u_r^2 + u_i^2 - v^2)*τ*u_r,
                          D(u_i) ~  ω * u_r + (u_r^2 + u_i^2 - v^2)*τ*u_i],
                         [ω, v], [u_i, u_r])

# #### Droop control

@parameters K u_ref x_ref x(t)
@variables u(t)

droop_control = IOBlock([
    u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
    ], [x], [u])

#=
After defining the blueprints we can create the Blocks of the system
based on them:
=#

p_filter = IOBlock(lpf, name = :p_filter)
q_filter = IOBlock(lpf, name = :q_filter)
p_droop = IOBlock(droop_control, name = :p_droop)
q_droop = IOBlock(droop_control, name = :q_droop)
v_source = IOBlock(voltage_source, name = :v_source)

nothing #hide
#=
We can put the blocks together to form this system:

```
     +----------+  +-------------+
 P --| p_filter |--|   p_droop   |   +----------+
     |    τ     |  | u_ref, xref |---|          |
     +----------+  +-------------+  ω| v_source |-- u_r
                                     |    τ     |
     +----------+  +-------------+  v|          |-- u_i
 Q --| q_filter |--|   q_droop   |---|          |
     |    τ     |  | u_ref, xref |   +----------+
     +----------+  +-------------+

```
=#
gfi = IOSystem([p_filter.filtered => p_droop.x,
                q_filter.filtered => q_droop.x,
                p_droop.u => v_source.ω,
                q_droop.u => v_source.v],
               [p_filter, q_filter, p_droop, q_droop, v_source],
               name = :GridForming,
               namespace_map = [p_filter.input => :P, q_filter.input => :Q, q_droop.u => :v],
               outputs = [v_source.u_i, v_source.u_r])

#=
There is still one problem: our connected system has the inputs `P` and `Q`. In order to construct
IONodes the IOSystems have to be systems which convert some complex current to complex voltage:
```
       +---------+
i_r -->|   ???   |--> u_r
i_i -->|         |--> u_i
       +---------+
```

We can achieve this by defining a nother block which converts `(i, u) ↦ (P, Q)`

```
               +----------+  +-------------+
     +-----+   | p_filter |--|   p_droop   |   +----------+
i_r--|     |-P-|    τ     |  | u_ref, xref |---|          |
i_i--| pow |   +----------+  +-------------+  ω| v_source |---+--u_r
     |     |                                   |    τ     |   |
 +---|     |   +----------+  +-------------+  v|          |-+-|--u_i
 | +-|     |-Q-| q_filter |--|   q_droop   |---|          | | |
 | | +-----+   |    τ     |  | u_ref, xref |   +----------+ | |
 | |           +----------+  +-------------+                | |
 | +--------------------------------------------------------+ |
 +-----------------------------------------------------------+

```
=#
@parameters u_i(t) u_r(t) i_i(t) i_r(t)
@variables P(t) Q(t)
pow = IOBlock([P ~ u_r*i_r + u_i*i_i,
                 Q ~ u_i*i_r - u_r*i_i],
                [u_i, u_r, i_i, i_r], [P, Q], name=:pow)

gfi2 = IOSystem([gfi.u_i => pow.u_i,
                 gfi.u_r => pow.u_r,
                 pow.P => gfi.P,
                 pow.Q => gfi.Q],
                [gfi, pow], outputs=[gfi.u_i, gfi.u_r], name=:gfi_i_to_u)
nothing #hide

# and the system can be reduced to a single IOBlock
connected = connect_system(gfi2)
nothing #hide

#=
The following part should live inside `PowerDynamics.jl` at some point. It defines an
`IONode` which will be constructed from an `IOBlock` as a new `AbstractNode` type.
=#
using PowerDynamics
using NetworkDynamics: ODEVertex
struct IONode{T} <: AbstractNode
    block::IOBlock
    generated::T
    parameters::Vector{Float64}
end

function IONode(blk::IOBlock, parameters::Dict)
    ## BlockSpec: blk must by of type (i_r, i_i) ↦ (u_r, u_i)
    spec = BlockSpec([:i_r, :i_i], [:u_r, :u_i])
    @assert spec(blk) "Block has to follow PowerDynamics i/o conventions!"
    ## TODO check parameters
    gen = generate_io_function(blk,
                               f_states=[blk.u_r, blk.u_i],
                               f_inputs=[blk.i_r, blk.i_i],
                               f_params=keys(parameters));
    IONode(blk, gen, collect(values(parameters)))
end

import PowerDynamics.construct_vertex
function construct_vertex(ion::IONode)
    gen = ion.generated
    function rhs!(dx, x, e_s, e_d, _p, t)
        i = total_current(e_s, e_d)
        gen.f_ip(dx, x, [real(i), imag(i)], ion.parameters, t)
    end
    ODEVertex(f! = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = Symbol.(gen.states))
end

import PowerDynamics.symbolsof
PowerDynamics.symbolsof(ionode::IONode) = Symbol.(ionode.generated.states)
PowerDynamics.dimension(ionode::IONode) = length(ionode.generated.states)


#=
Let's try it out:
The parameters have to be provided as a Dict. This dict should contain all internal parameters of the block.
The list of parameters can be retrieved from the block:
=#
connected.iparams

#=
The parameter dict can be defined in several ways, I think the first option is the most convenient
because it is a pure `Symbol` which does not depend on any variables in the scope. The namespace separator
can be typed as `\_+<tab>`.
```
:q_filter₊τ => 1.0
q_filter.τ => 1.0
:gfi_i_to_u₊q_filter₊τ => 1.0
connected.q_filter₊τ => 1.0
```
=#
para = Dict(:p_filter₊τ => 1.0, # super legit parameter choice
            :q_filter₊τ => 1.0,
            :v_source₊τ => 1.0,
            :p_droop₊u_ref => 1.0,
            :q_droop₊u_ref => 1.0,
            :p_droop₊x_ref => 1.0,
            :q_droop₊x_ref => 1.0,
            :q_droop₊K => 1.0,
            :p_droop₊K => 1.0)
n = IONode(connected, para)
nothing #hide

# Now let's test whether this works in PD..
using PowerDynamics: SlackAlgebraic, PiModelLine,PowerGrid,find_operationpoint
using OrderedCollections: OrderedDict

buses=OrderedDict(
    "bus1"=> IONode(connected, para),
    "bus2"=> SlackAlgebraic(U=1))

branches=OrderedDict(
    "branch1"=> PiModelLine(from= "bus1", to = "bus2",y=4.999131600798035-1im*15.263086523179553, y_shunt_km=0.0528/2, y_shunt_mk=0.0528/2))

powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid)
timespan= (0.0,5.)
nothing #hide

# simulating a voltage perturbation at node 1
fault1 = ChangeInitialConditions(node="bus1", var=:u_r, f=Inc(0.2))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
using Plots
plot(solution1.dqsol)

# well, it does something.
