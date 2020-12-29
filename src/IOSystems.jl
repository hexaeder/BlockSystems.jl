module IOSystems

using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: Parameter, ODESystem, Differential
using ModelingToolkit: rename, getname, renamespace, namespace_equations
using ModelingToolkit: equation_dependencies, asgraph, variable_dependencies, eqeq_dependencies, varvar_dependencies
using SymbolicUtils: Symbolic, to_symbolic
using LightGraphs

export AbstractIOSystem, IOBlock, IOSystem, connect_system

@doc raw"""
    AbstractIOSystem

Basic type for IOSystems. Such systems contains a ODEProblem of the form
```math
\begin{aligned}
\dot \mathbf x(t) &= f( \mathbf i(t) )\\
\mathbf o(t) &= g(\mathbf x(t), \mathbf i(t))
\end{aligned}
```
where the states can be separated into `states = inputs ∪ istates ∪ outputs` (`istates` is short for `internalstates`).

Each `AIOS` has a `name` field generating a namespace. The symbols of the
"""
abstract type AbstractIOSystem end

function Base.getproperty(sys::AbstractIOSystem, name::Symbol)
    if name ∈ fieldnames(typeof(sys))
        return getfield(sys, name)
    end
    syms = vcat(sys.inputs, sys.iparams, sys.istates, sys.outputs)

    i = findfirst(x -> getname(x) == name, syms)
    if i !== nothing
        return rename(syms[i], renamespace(sys.name, name))
    end
    throw(error("Variable $name does not exist"))
end

function namespace_property(ios::AbstractIOSystem, prop::Symbol)
    [rename(sym, renamespace(ios.name, getname(sym)))
     for sym in getproperty(ios, prop)]
end
namespace_inputs(ios::AbstractIOSystem)  = namespace_property(ios, :inputs)
namespace_iparams(ios::AbstractIOSystem) = namespace_property(ios, :iparams)
namespace_istates(ios::AbstractIOSystem) = namespace_property(ios, :istates)
namespace_outputs(ios::AbstractIOSystem) = namespace_property(ios, :outputs)

isunique(collection) = length(collection) == length(unique(collection))

"""
$(TYPEDEF)

A basic IOSystem which consists of a single ODESystems.
$(TYPEDFIELDS)
"""
struct IOBlock <: AbstractIOSystem
    name::Symbol
    inputs::Vector{Symbolic}
    iparams::Vector{Symbolic}
    istates::Vector{Symbolic}
    outputs::Vector{Symbolic}
    system::ODESystem
end

function IOBlock(eqs::AbstractVector{<:Equation}, inputs, outputs; name = gensym(:IOBlock))
    os = ODESystem(eqs, name = name)

    # TODO: check constrsaints to system of equations
    # TODO: remove assertions in favour of errors
    @assert Set(inputs) ⊆ Set(parameters(os)) "inputs musst be parameters"
    @assert Set(outputs) ⊆ Set(states(os)) "outputs musst be variables"

    inputs = parameters(os) ∩ inputs # gets the inputs as `tern` type
    outputs = os.states ∩ outputs # gets the outputs as `tern` type
    istates = setdiff(os.states, outputs)
    iparams = setdiff(parameters(os), inputs)

    @assert Set(os.states) == Set(outputs ∪ istates)
    @assert Set(parameters(os)) == Set(inputs ∪ iparams)

    IOBlock(name, inputs, iparams, istates, outputs, os)
end

"""
$(TYPEDEF)

A composit IOSystem which consists of multipe `AbstractIOSystem`(@ref) wich are connected via
a vector of namespaced pairs (`subsys1.out1 => subsys2.in1`).
The inputs and outputs of the composit system are declared in the namespace of the new composit
system but point to namespaced variables of the subsystems.

$(TYPEDFIELDS)
"""
struct IOSystem <: AbstractIOSystem
    name::Symbol
    inputs::Vector{Symbolic}
    iparams::Vector{Symbolic}
    istates::Vector{Symbolic}
    outputs::Vector{Symbolic}
    connections::Dict{Symbolic, Symbolic}
    inputs_map::Dict{Symbolic, Symbolic}
    iparams_map::Dict{Symbolic, Symbolic}
    istates_map::Dict{Symbolic, Symbolic}
    outputs_map::Dict{Symbolic, Symbolic}
    systems::Vector{AbstractIOSystem}
end

function IOSystem(cons,
                  io_systems::Vector{<:AbstractIOSystem};
                  inputs_map = nothing,
                  iparams_map = nothing,
                  istates_map = nothing,
                  outputs_map = nothing,
                  name = gensym(:IOSystem),
                  )
    namespaces = [sys.name for sys in io_systems]
    @assert namespaces == unique(namespaces) "Namespace collision in subsystems!"

    @assert isunique(first.(cons)) "Multiple connections to same input!"
    namespaced_inputs = vcat([namespace_inputs(sys)
                              for sys in io_systems]...)
    namespaced_iparams = vcat([namespace_iparams(sys)
                               for sys in io_systems]...)
    namespaced_istates = vcat([namespace_istates(sys)
                               for sys in io_systems]...)
    namespaced_outputs = vcat([namespace_outputs(sys)
                               for sys in io_systems]...)

    @assert Set(first.(cons)) ⊆ Set(namespaced_inputs) "First argument in connection needs to be input of subsystem."
    @assert Set(last.(cons)) ⊆ Set(namespaced_outputs) "Second argument in connection needs to be output of subsystem."

    # namespace promotion for inputs
    if inputs_map === nothing
        inputs_map = create_namespace_map(io_systems, :inputs, skip = first.(cons))
    else
        open_inputs = setdiff(namespaced_inputs, first.(cons))
        inputs_map = fix_map_types(inputs_map)
        @assert keys(inputs_map) ⊆ Set(open_inputs) "inputs_map !⊆ open_inputs"
        @assert isunique(values(inputs_map)) "rhs of inputs_map not unique"
        # autonamespace all the other inputs
        skip = Set(keys(inputs_map)) ∪ first.(cons)
        remaining = create_namespace_map(io_systems, :inputs, skip = skip)
        inputs_map = merge(inputs_map, remaining)
    end

    # namespace promotion for iparams
    if iparams_map === nothing
        iparams_map = create_namespace_map(io_systems, :iparams)
    else
        iparams_map = fix_map_types(iparams_map)
        @assert keys(iparams_map) ⊆ Set(namespaced_iparams) "iparams_map !⊆ iparams"
        @assert isunique(values(iparams_map)) "rhs of iparams_map not unique"
        # autonamespace all the other inputs
        skip = Set(keys(iparams_map))
        remaining = create_namespace_map(io_systems, :iparams, skip = skip)
        iparams_map = merge(iparams_map, remaining)
    end

    # namespace promotion for istates
    if istates_map === nothing
        istates_map = create_namespace_map(io_systems, :istates)
    else
        istates_map = fix_map_types(istates_map)
        @assert keys(istates_map) ⊆ Set(namespaced_istates) "istates_map !⊆ istates"
        @assert isunique(values(istates_map)) "rhs of istates_map not unique"
        # autonamespace all the other inputs
        skip = Set(keys(istates_map))
        remaining = create_namespace_map(io_systems, :istates, skip = skip)
        istates_map = merge(istates_map, remaining)
    end

    # namespace promotion for outputs
    if outputs_map === nothing
        outputs_map = create_namespace_map(io_systems, :outputs)
    else
        outputs_map = fix_map_types(outputs_map)
        @assert keys(outputs_map) ⊆ Set(namespaced_outputs) "outputs_map !⊆ outputs"
        @assert isunique(values(outputs_map)) "rhs of outputs_map not unique"
    end

    # the automatic namespace promotion is not aware of the names used
    # in the *_maps given by the user. Therfore a namespace clash might
    # happen. This won't be the case for outputs since there is no remaining
    # part which is promoted automatically.
    @assert isunique(values(inputs_map)) "namespace promotion of inputs clashed with manually given inputs_map"
    @assert isunique(values(iparams_map)) "namespace promotion of iparams clashed with manually given inputs_map"
    @assert isunique(values(istates_map)) "namespace promotion of istates clashed with manually given inputs_map"

    cons = Dict(cons)

    IOSystem(
        name,
        collect(values(inputs_map)),
        collect(values(iparams_map)),
        collect(values(istates_map)),
        collect(values(outputs_map)),
        cons,
        inputs_map,
        iparams_map,
        istates_map,
        outputs_map,
        io_systems,
      )
end

"""
    fix_map_types(map)

Creates Dict from map. Changes `Num`-types to `Symbolic`-types in the
user provided maps. Map can be Dict or  Array of Pairs
"""
function fix_map_types(map)
    dict = Dict(map)
    newdict = Dict{Symbolic, Symbolic}()
    for k in keys(dict)
        newkey = to_symbolic(k)
        newdict[newkey] = to_symbolic(dict[k])
    end
    newdict
end

"""
    create_namespace_map(io_systems, property; skip =[])

Creates dictionary from namespaced=>promoted symbols of property.

property ∈ {:inputs, :iparams, :istates, :outputs}

Tries to get rid of the namespace wherever possible. Namespaced symbols
which should not be included in the resulting map can be specified using
skip parameter.
"""
function create_namespace_map(io_systems, property; skip = [])
    all = vcat((getproperty(ios, property) for ios in io_systems)...)
    all_spaced = vcat((namespace_property(ios, property) for ios in io_systems)...)

    # namespaced => non-namespaced pairs
    pairs = [Pair(t...) for t in zip(all_spaced, all)]
    # filter out values which should be skipped (i.e. connected inputs)
    filter!(x -> first(x) ∉ Set(skip), pairs)

    # find duplicates and add namepsaces
    duplicates = Set()
    for i in 1:length(pairs)
        if last(pairs[i]) ∈ Set(last.(pairs[i+1:end]))
            push!(duplicates, last(pairs[i]))
        end
    end
    # no namespace promotion for duplicates
    for i in 1:length(pairs)
        if pairs[i].second ∈ duplicates
            pairs[i] = pairs[i].first => pairs[i].first
            @warn "Could not promote $(pairs[i].first) to system namespace."
        end
    end

    return Dict(pairs)
end

function connect_system(ios::IOSystem)
    # recursive connect all subsystems
    for (i, subsys) in enumerate(ios.systems)
        if subsys isa IOSystem
            ios.systems[i] = connect_system(subsys)
        end
    end
    eqs = vcat([namespace_equations(iob.system) for iob in ios.systems]...)

    # get rid of closed inputs by substituting output states
    connections = ios.connections
    for (i, eq) in enumerate(eqs)
        eqs[i] = eq.lhs ~ substitute(eq.rhs, connections)
    end

    # TODO: poove assumtions
    # - each state is represented by one lhs
    # - lhs only first order or algebraic
    # - no self loop in algeraic

    # reduce algebraic states of the system
    eqs = reduce_algebraic_states(eqs, skip = keys(ios.outputs_map))

    # TODO: possibly get red of unused states
    # (internal variables which are not used for the outputs)
    # hint: attracting_components(dependency graph)
    # gets partions whitout leaving edges which means,
    # this subset of equations is NOT used by the others

    # apply the namespace transformations
    namespace_promotion = merge(ios.inputs_map, ios.iparams_map, ios.istates_map, ios.outputs_map)

    for (i, eq) in enumerate(eqs)
        eqs[i] = eqsubstitute(eq, namespace_promotion)
    end

    # eqs
    IOBlock(eqs, ios.inputs, ios.outputs, name=ios.name)
end

function reduce_algebraic_states(eqs::Vector{Equation}; skip=[])
    eqs = deepcopy(eqs)
    sys = ODESystem(eqs) # will be used for the dependency graph
    eqs = sys.eqs # the ODESystem might reorder the equations

    # only consider states for reduction which are explicit algebraic and not in skip
    condition = eq -> is_explicit_algebraic(eq) && !(Set(get_variables(eq.lhs)) ⊆ Set(skip))
    algebraic_idx = findall(condition, eqs)

    # generate dependency graph
    graph = eqeq_dependencies(asgraph(sys), variable_dependencies(sys))

    # algebraic states which do no have cyclic dependencies between them can be reduced
    cycles = simplecycles(graph)
    reducable = pairwise_cycle_free(algebraic_idx, cycles)
    rules = [eq.lhs => eq.rhs for eq in eqs[reducable]]
    # subsitute all the equations, remove substituted
    for (i, eq) in enumerate(eqs)
        eqs[i] = eq.lhs ~ substitute(eq.rhs, rules)
    end
    deleteat!(eqs, reducable)
    return eqs
end

function pairwise_cycle_free(idx, cycles)
    pairs = collect(Iterators.product(idx, idx))
    # check for each pair whether found in cycle
    incycle(p) = p[1]==p[2] ? false : any(map(c -> p ⊆ c, cycles))
    A = incycle.(pairs)
    if !any(A)
        return idx
    else
        # find the variable with the highes number of loops
        cycle_count = reduce(+, A, dims=1)[:]
        worst_idx = findmax(cycle_count)[2]
        deleteat!(idx, worst_idx)
        return pairwise_cycle_free(idx, cycles)
    end
end

function is_explicit_algebraic(eq::Equation)
    if eq.lhs isa Sym || eq.lhs isa Term && !(eq.lhs.op isa Differential)
        vars = get_variables(eq.lhs)
        @assert length(vars) == 1
        return vars[1] ∉ Set(get_variables(eq.rhs))
    else
        return false
    end

end

eqsubstitute(eq::Equation, rules) = substitute(eq.lhs, rules) ~ substitute(eq.rhs, rules)

"""
The substitution only works if the Terms have the same eltype.
"""
fixtermtype(sym::Sym) = sym
fixtermtype(term::Term{Real}) = term
fixtermtype(term::Term{Parameter{Real}}) = Term{Real}(term.op, term.args)
fixtermtype(pair::Pair{Symbolic, Symbolic}) = fixtermtype(pair.first)=>fixtermtype(pair.second)
function fixtermtype(dict::Dict{Symbolic, Symbolic})
    new = Dict{Symbolic, Symbolic}()
    for pair in dict
        push!(new, fixtermtype(pair))
    end
    new
end

# TODO: during namespace_equations the Term{Parameter{Real}} -> Term{Real} while the Sym{Parameter{Real}} stay at it is => PR in ModelingToolkit


end
