var documenterSearchIndex = {"docs":
[{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"EditURL = \"https://github.com/hexaeder/IOSystems_prototype/blob/master/examples/spacecraft.jl\"","category":"page"},{"location":"generated/spacecraft/#Control-System:-Spaceship","page":"Spacecraft","title":"Control System: Spaceship","text":"","category":"section"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"As an example we want to build an IOSystem controlling the altitude of a spacecraft. The spacecraft has mass m and can be controlled with thrusters which apply the force F(t) to the spacecraft. The altitude x(t)","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"dot v(t) = frac F(t) m\ndot x(t) =  v(t)","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"in our model this system has the input F(t), the internal state v(t) (vertical velocity) and the output x(t).","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"       *------------*\nF(t) --| spacecraft |-- x(t)\n       | m, v(t)    |\n       *------------*\n","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"using IOSystems\nusing ModelingToolkit\n@parameters t M F(t)\n@variables x(t) v(t)\n@derivatives D'~t\n\nspacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v], # define the equation\n                     [F], # inputs of the system\n                     [x], # outputs of the system\n                     name = :spacecraft)\n@show spacecraft.inputs\n@show spacecraft.outputs\n@show spacecraft.istates\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"We want to model a controller which takes a desired altitude as an input parameter and outputs the force for thrusters.","category":"page"},{"location":"generated/spacecraft/#Simple-proportional-controller","page":"Spacecraft","title":"Simple proportional controller","text":"","category":"section"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"A proportional controller takes an input i and calculates the output proportional to the input.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":" o(t) = Kcdot i(t)","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"       *--------*\ni(t) --| prop K |-- o(t)\n       *--------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"@parameters K i(t)\n@variables o(t)\n\nprop = IOBlock([o ~ K*i], [i], [o], name = :prop)\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"In order to make this useful as an controller, the input has to be the difference between the reference and the system variable (negative feedback). We can model this as an IOSystem where","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Δ = p - m","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"       *--------*\np(t) --|  diff  |-- Δ(t)\nm(t) --|        |\n       *--------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"@parameters p(t) m(t)\n@variables Δ(t)\ndiff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Now we can connect both of the defined models to create an proportional controller","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"             *----------------------------------------*\n             | propc                                  |\n             |         *--------*   *--------*        |\n  target(t)--|--p(t) --|  diff  |---| prop K |--o(t)--|--o(t)\nfeedback(t)--|--m(t) --|        |   *--------*        |\n             |         *--------*                     |\n             *----------------------------------------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"If we don't provide additional information the system will try to promote all of the enclosed variables to the new syste mwide namespace.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"prop_c = IOSystem([prop.i => diff.Δ], # connect input of prop to output of diff\n                  [diff, prop], # subsystems\n                  name=:propc)\n@info \"namespace mapping\" prop_c.inputs_map prop_c.istates_map prop_c.outputs_map","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"For finer control it is often prefered to give new names manually:","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"@parameters target(t) feedback(t)\nprop_c = IOSystem([prop.i => diff.Δ], [diff, prop],\n                  outputs_map = [prop.o => o],\n                  inputs_map = [diff.p => target, diff.m => feedback],\n                  name=:propc)\n@info \"namespace mapping\" prop_c.inputs_map prop_c.istates_map prop_c.outputs_map","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Right now, the created object is a container for the two included systems. However, it is possible to transform the object into a new IOBlock by calling the connect_system function. The resulting is equivalent to","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"             *----------------------------------*\n  target(t)--| prop_c_block                     |--o(t)\nfeedback(t)--| o(t)=K*(target(t) - feedback(t)) |\n             *----------------------------------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"prop_c_block = connect_system(prop_c)\nnothing #hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Now we can hook our spaceship to this controller. It does not matter whether we use the connected IOBlock version prop_c or the IOSystem version prop_c_block. We want to build the connected system","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"           *--------------------------------------------*\n           |  control system                            |\n           |       *--------*   *------------*          |\ntarget(t)--|-------| prop_c |---| spacecraft |-x(t)--*--|--altitude(t)\n           |  *-fb-|        |   | m, v(t)    |       |  |\n           |  |    *--------*   *------------*       |  |\n           |  *--------------------------------------*  |\n           *--------------------------------------------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"@variables altitude(t)\nspace_controller = IOSystem([spacecraft.F => prop_c.o, prop_c.feedback => spacecraft.x],\n                            [prop_c, spacecraft],\n                            outputs_map = [spacecraft.x => altitude])\n# we want to reduce the space_controller to a block\nspace_controller = connect_system(space_controller)\n@info \"Variables of space_controller\" space_controller.inputs space_controller.outputs space_controller.istates space_controller.iparams space_controller.system.eqs","category":"page"},{"location":"generated/spacecraft/#Simulate-System","page":"Spacecraft","title":"Simulate System","text":"","category":"section"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"In order to simulate the system we can have to build the Julia functions.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"gen = generate_io_function(space_controller)\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"By doing so we get access to a named tuple with the fileds","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"gen.f_ip in-place function\ngen.f_oop out-of-place function\ngen.massm mass matrix of the system\ngen.states symbols of states (in order)\ngen.inputs symbols of inputs (in order)\ngen.params symbols of parameters (in order)","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"The functions have the form f_ip(du, u, inputs, params, t) where u are all the states (outputs stacked on top of internal states) and t is the independent variable of the system. The order of the inputs and states can be controlled.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"gen = generate_io_function(space_controller, first_states=[altitude])\n@info \"Generated function\" gen.massm gen.states gen.inputs gen.params\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Well, let's see how our model is doing.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"using Plots\nusing DifferentialEquations\ntargetfun(t) = t>1.0 ? 1.0 : 0\nodefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)\np = [0.5, 1.0] # K, m\nu0 = [0.0, 0.0] # altitude, v\ntspan = (0.0, 30.0)\nprob = ODEProblem(odefun, u0, tspan, p)\nsol = solve(prob)\nplot(t->sol(t)[1],tspan..., label=\"altitude\", title=\"proportional control\")\nplot!(t->targetfun(t),tspan..., label=\"target\")","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"Well who could have thought, proportional control looks like an harmonic oscillator 🤷‍♂️","category":"page"},{"location":"generated/spacecraft/#Defining-a-better-controller","page":"Spacecraft","title":"Defining a better controller","text":"","category":"section"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"We might just add a damping term (a force proportional to the velocity of the spaceship). If it works for a harmonic oscillator, it should work for our spaceship.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"           *--------------------------------------------------------*\n           |  control system                                        |\n           |  *--------------------------------------------------*  |\n           |  |    *--------*     *---*                          |  |\n           |  *-v--| prop_v |-(-)-| d |     *------------*       |  |\n           |       *--------*     | i |--F--| spacecraft |-v(t)--*  |\n           |       *--------*     | f |     | m          |-x(t)--*--|--altitude(t)\ntarget(t)--|-------| prop_c |-(+)-| f |     *------------*       |  |\n           |  *-fb-|        |     *---*                          |  |\n           |  |    *--------*                                    |  |\n           |  *--------------------------------------------------*  |\n           *--------------------------------------------------------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"In order to do so we have to slightly redefine the spaceship system: now the velocity v(t) is also an output and not and internal state.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v],\n                     [F],\n                     [x,v],\n                     name = :spacecraft)","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"One can define new blocks based on previously defined blocks.","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"prop_v = IOBlock(prop, name=:prop_v)\nfdiff = IOBlock(diff, name=:fdiff)\n\nspace_controller = IOSystem([prop_v.i => spacecraft.v,\n                             prop_c.feedback => spacecraft.x,\n                             fdiff.p => prop_c.o,\n                             fdiff.m => prop_v.o,\n                             spacecraft.F => fdiff.Δ],\n                            [prop_v, prop_c, fdiff, spacecraft],\n                            outputs_map = [spacecraft.x => altitude])\n\nspace_controller = connect_system(space_controller)\ngen = generate_io_function(space_controller, first_states=[altitude])\n\nodefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)\np = [1.0, 1.0, 0.5] # propc, M, propv\nu0 = [0.0, 0.0] # altitude, v\ntspan = (0.0, 30.0)\nprob = ODEProblem(odefun, u0, tspan, p)\nsol = solve(prob)\nplot(sol, vars=(0,1), label=\"altitude\", title=\"better control\")\nplot!(t->targetfun(t),tspan..., label=\"target\")","category":"page"},{"location":"generated/spacecraft/#Defining-an-PT1-controller","page":"Spacecraft","title":"Defining an PT1 controller","text":"","category":"section"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"             *-----------------------------------------------*\n             | pi_c                                   *---*  |\n             |                           *------------|   |  |\n             |        *----*  *--------* | *-------*  |sum|--|--o(t)\n  target(t)--|--p(t)--|diff|--| prop K |-*-| int T |--|   |  |\nfeedback(t)--|--m(t)--|    |  *--------*   *-------*  *---*  |\n             |        *----*                                 |\n             *-----------------------------------------------*","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"@parameters T a1(t) a2(t)\n@variables Σ(t) altitude(t)\nint = IOBlock([D(o) ~ 1/T * i - o], [i], [o], name=:int)\nadder = IOBlock([Σ ~ a1 + a2], [a1, a2], [Σ], name=:add)\n\npi_c = IOSystem([prop.i => diff.Δ,\n                 int.i => prop.o,\n                 adder.a1 => prop.o,\n                 adder.a2 => int.o],\n                [diff, prop, int, adder],\n                outputs_map = [adder.Σ => o],\n                inputs_map = [diff.p => target, diff.m => feedback],\n                name=:pi_c)\nnothing # hide","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"as before we can close the loop and build the control circuit","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"space_controller = IOSystem([spacecraft.F => pi_c.o, pi_c.feedback => spacecraft.x],\n                            [pi_c, spacecraft],\n                            outputs_map = [spacecraft.x => altitude])\nspace_controller = connect_system(space_controller, verbose=false)\n@info \"Variables of space_controller\" space_controller.inputs space_controller.outputs space_controller.istates space_controller.iparams space_controller.system.eqs","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"and we can simulate and plot the system","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"gen = generate_io_function(space_controller, first_states=[altitude])\n@info \"order of parameters\" gen.params\n\nodefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)\np = [0.5, -1.5, 1.0] # K, T, m\nu0 = [0.0, 0.0, 0.0] # altitude, int.o, v\ntspan = (0.0, 50.0)\nprob = ODEProblem(odefun, u0, tspan, p)\nsol = solve(prob)\nplot(sol, vars=(0,[ 1,2 ]), label=[\"altitude\" \"integrator\"], title=\"PI controller\")\nplot!(t->targetfun(t),tspan..., label=\"target\")","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"thank you for flying with us :) TODO: FIX ordering of parameters, the results change from doc build to doc build","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"","category":"page"},{"location":"generated/spacecraft/","page":"Spacecraft","title":"Spacecraft","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"EditURL = \"https://github.com/hexaeder/IOSystems_prototype/blob/master/examples/kuramoto_network.jl\"","category":"page"},{"location":"generated/kuramoto_network/#Integration-with-NetworkDynamics.jl","page":"Kuramoto Network","title":"Integration with NetworkDynamics.jl","text":"","category":"section"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"In this example we model a Network based on the Kuroamoto model using the NetworkDynamics.jl package.","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"using NetworkDynamics\nusing LightGraphs\nusing IOSystems\nusing ModelingToolkit\nusing DifferentialEquations\nusing Plots","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"In the Kuramoto model each vertex i has it's own intrinsic angular frequency omega_i. The angle of frequency is given by","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"dotphi = omega_i + sum e","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"where sum e sums over all the connected edges.","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"@parameters t ω edgesum(t)\n@variables ϕ(t)\n@derivatives D'~t\n\nkmvert = IOBlock([D(ϕ) ~ ω + edgesum],\n                 [edgesum],\n                 [ϕ],\n                 name=:kmvert)\n\nfunction generate_ode_vertex(iob::IOBlock, aggregator; kwargs...)\n    gen = generate_io_function(iob; kwargs...)\n    f = (du, u, edges, p, t) -> gen.f_ip(du, u, aggregator(edges), p, t)\n    vars = ModelingToolkit.tosymbol.(gen.states)\n    ODEVertex(f! = f, dim = length(vars), sym = vars)\nend\n\n# allocation free oop aggregator. might be more difficult for more-dimensional systems\n# unfortunately there are no vector variables in MDK and we can't model the aggregators\n# as IOSystems.\nfunction sumedges(edges)\n    r = 0.0\n    for e in edges\n        r += e[1]\n    end\n    return r\nend\n\nvertex = generate_ode_vertex(kmvert, sumedges)\nnothing #hide","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"The edges are defined as","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"e_ij = K cdot mathrmsin(phi_i - phi_j)","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"For the (static) edges we can't used IOSystem function building in the moment. Well, we could but than we'd get an ODEEdge with zero-massmatrix.","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"@parameters v₁(t) v₂(t) K\nedge_ip = build_function([K*sin(v₁[1] - v₂[1])], v₁, v₂, K, t, expression=Val{false})[2]\nedge = StaticEdge(f! = edge_ip, dim = 1, coupling=:antisymmetric)\n\nnothing #hide","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"Lets runt the same example as in NetworkDynamics.jl","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"N = 8\ng = watts_strogatz(N,2,0) # ring network\nnd = network_dynamics(vertex, edge, g)\n\nω = collect(1:N)./N\nω .-= sum(ω)/N\nK = 3.\np = (ω, K); # p[1] vertex parameters, p[2] edge parameters\n\nx0 = collect(1:N)./N\nx0 .-= sum(x0)./N\ntspan = (0., 10.)\nprob = ODEProblem(nd, x0, tspan, p)\nsol = solve(prob)\nplot(sol, ylabel=\"ϕ\")","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"","category":"page"},{"location":"generated/kuramoto_network/","page":"Kuramoto Network","title":"Kuramoto Network","text":"This page was generated using Literate.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = IOSystems","category":"page"},{"location":"#IOSystems","page":"Home","title":"IOSystems","text":"","category":"section"},{"location":"#Basics","page":"Home","title":"Basics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"An input-output-system is characterized by a set of equations. These equations can be either first order ODEs or explicit algebraic equations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\ndot mathbf x(t) = f(mathbf x(t) mathbf y(t) mathbf i(t) p)\nmathbf y(t) = g(mathbf x(t) mathbf y(t) mathbf i(t) p)\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"such system contains of ","category":"page"},{"location":"","page":"Home","title":"Home","text":"states (x and y) and\nparameters (i, p).","category":"page"},{"location":"","page":"Home","title":"Home","text":"States are determined by the given equations (i.e. the equations describe how the states change). Parameters are externally given. For IO systems we define subgroups","category":"page"},{"location":"","page":"Home","title":"Home","text":"states\ninternal states (istates) which are meant for internal use\noutput states (outputs) which might be used as inputs for other systems\nparameters\ninternal parameters (iparams) which are typically constant and\ninputs (inputs) which can be connected to the outputs of other systems.","category":"page"},{"location":"#Types","page":"Home","title":"Types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The base type for the IOSystems is AbstractIOSystem with the has two concrete implementations. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"AbstractIOSystem","category":"page"},{"location":"#IOSystems.AbstractIOSystem","page":"Home","title":"IOSystems.AbstractIOSystem","text":"abstract supertype for IOBlock and IOSystem.\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"An IOBlock consists of a set of equations, a set of inputs and outputs and a name.","category":"page"},{"location":"","page":"Home","title":"Home","text":"IOBlock\nIOBlock(eqs::AbstractVector{<:Equation}, inputs, outputs; name = gensym(:IOBlock))","category":"page"},{"location":"#IOSystems.IOBlock","page":"Home","title":"IOSystems.IOBlock","text":"struct IOBlock <: AbstractIOSystem\n\nA basic IOSystem which consists of a single ODESystem.\n\nname::Symbol\ninputs::Array{SymbolicUtils.Symbolic,1}\niparams::Array{SymbolicUtils.Symbolic,1}\nistates::Array{SymbolicUtils.Symbolic,1}\noutputs::Array{SymbolicUtils.Symbolic,1}\nsystem::ODESystem\n\n\n\n\n\n","category":"type"},{"location":"#IOSystems.IOBlock-Tuple{AbstractArray{var\"#s23\",1} where var\"#s23\"<:Equation,Any,Any}","page":"Home","title":"IOSystems.IOBlock","text":"IOBlock(eqs, inputs, outputs; name)\n\n\nConstruct a new IOBlock for the given arguments.\n\nusing IOSystems, ModelingToolkit\n@parameters t i(t)\n@variables x(t) o(t)\n@derivatives D'~t\n\niob = IOBlock([D(x) ~ i, o ~ x], [i], [o], name=:iob)\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"An IOSystem consists of multiple AbstractIOSystems and the connections between them.","category":"page"},{"location":"","page":"Home","title":"Home","text":"IOSystem\nIOSystem(cons, io_systems::Vector{<:AbstractIOSystem}; inputs_map = nothing, iparams_map = nothing, istates_map = nothing, outputs_map = nothing, name = gensym(:IOSystem), autopromote = true)","category":"page"},{"location":"#IOSystems.IOSystem","page":"Home","title":"IOSystems.IOSystem","text":"struct IOSystem <: AbstractIOSystem\n\nA composite IOSystem which consists of multiple AbstractIOSystem which are connected via a vector of namespaced pairs (subsys1.out1 => subsys2.in1).\n\nAn IOSystem contains maps how to promote the namespaced variables of the subystem to the new scope    subsys1₊x(t) => x(t)    subsys1₊y(t) => subsys1₊y(t)    subsys2₊y(t) => subsys2₊y(t)\n\nname::Symbol\ninputs::Array{SymbolicUtils.Symbolic,1}\niparams::Array{SymbolicUtils.Symbolic,1}\nistates::Array{SymbolicUtils.Symbolic,1}\noutputs::Array{SymbolicUtils.Symbolic,1}\nconnections::Dict{SymbolicUtils.Symbolic,SymbolicUtils.Symbolic}\ninputs_map::Dict{SymbolicUtils.Symbolic,SymbolicUtils.Symbolic}\niparams_map::Dict{SymbolicUtils.Symbolic,SymbolicUtils.Symbolic}\nistates_map::Dict{SymbolicUtils.Symbolic,SymbolicUtils.Symbolic}\noutputs_map::Dict{SymbolicUtils.Symbolic,SymbolicUtils.Symbolic}\nsystems::Array{AbstractIOSystem,1}\n\n\n\n\n\n","category":"type"},{"location":"#IOSystems.IOSystem-Tuple{Any,Array{var\"#s23\",1} where var\"#s23\"<:AbstractIOSystem}","page":"Home","title":"IOSystems.IOSystem","text":"IOSystem(cons, io_systems; inputs_map, iparams_map, istates_map, outputs_map, name, autopromote)\n\n\nConstruct a new IOSystem from various subsystems. Parameters\n\ncons: the connections in the form sub2.input => sub2.output\nio_systems: Vector of subsystems\ninputs_map, iparams_map, istates_map: Provide custom namespace promotion / renaming i.e. sub1.input => voltage. Variables without entry in the map will be promoted automatically.\noutputs_map: Provide custom namespace promotion for outputs. If this parameter is not given, all of the outputs (connected and unconnected) will be marked as outputs of the whole system. If the outputs map is given only the given outputs will be marked as outputs. All other outputs of subsystems will become internal states of the connected system (and might be optimized away in connect_system).\nname: namespace\nautopromote=true: enable/disable automatic promotion of variable names to system namespace\n\n\n\n\n\n","category":"method"},{"location":"#Transformations","page":"Home","title":"Transformations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Each IOSystem can be transformed into an IOBlock. At this step the actual manipulation of the equations happen.","category":"page"},{"location":"","page":"Home","title":"Home","text":"connect_system","category":"page"},{"location":"#IOSystems.connect_system","page":"Home","title":"IOSystems.connect_system","text":"connect_system(ios; verbose)\n\n\nRecursively transform IOSystems to IOBlocks.\n\nsubstitute inputs with connected outputs\ntry to eliminate equations for internal states which are not used to calculate the specified outputs of the system.\ntry to eliminate explicit algebraic equations (i.e. outputs of internal blocks) by substituting each occurence with their rhs. Explicit algebraic states which are marked as system outputs won't be reduced.\n\nParameters:\n\nios: system to connect\nverbose=false: toggle verbosity (show equations at different steps)\n\n\n\n\n\n","category":"function"},{"location":"#Function-building","page":"Home","title":"Function building","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"generate_io_function","category":"page"},{"location":"#IOSystems.generate_io_function","page":"Home","title":"IOSystems.generate_io_function","text":"generate_io_function(ios; first_states, first_inputs, simplify, expression, verbose)\n\n\nGenerate calleble functions for an AbstractIOSystem. An IOSystem will be transformed to an IOBlock first.\n\nParameters:\n\nios: the system to build the function\nfirst_states: define states=(istates ∪ outputs) which shoud appar first\nfirst_inputs: define inputs which shoud appar first\nsimplify=true: call simplify on the equations before building?\nexpression=Val{false}: toggle expression and callable function output\n\nReturns an named tuple with the fields\n\nf_ip in-place function f(dstates, states, inputs, params, iv)\nf_oop out-of-place function f(states, inputs, params, iv)=>dstates\nmassm mass matrix of the system\nstates symbols of states (in order)\ninputs symbols of inputs (in order)\nparams symbols of parameters (in order)\n\n\n\n\n\n","category":"function"}]
}
