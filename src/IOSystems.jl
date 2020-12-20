module IOSystems

using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: rename, getname, renamespace
using SymbolicUtils: Symbolic

export AbstractIOSystem, IOBlock, IOSystem

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

    # put asserts here
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

function IOSystem(
    cons,
    io_systems::Vector{<:AbstractIOSystem};
    inputs_map = nothing,
    outputs_map = nothing,
    name = gensym(:IOSystem),
)

    if inputs_map === nothing
        inputs_map =
            collect_and_resolve_namespace(io_systems, :inputs, skip = last.(cons))
        # for plug in last.(cons)
        #     delete!(inputs_map, findfirst(x->isequal(x, plug), inputs_map))
        # end
    end
    if outputs_map === nothing
        outputs_map = collect_and_resolve_namespace(io_systems, :outputs)
    end
    istates_map = collect_and_resolve_namespace(io_systems, :istates)

    cons = Dict(cons)

    IOSystem(
        name,
        collect(keys(inputs_map)),
        collect(keys(istates_map)),
        collect(keys(outputs_map)),
        cons,
        inputs_map,
        istates_map,
        outputs_map,
        io_systems,
    )
end

function create_namespace_map(io_systems, property; skip = [])
    all = vcat((getproperty(ios, property) for ios in io_systems)...)
    all_spaced = vcat((namespace_property(ios, property) for ios in io_systems)...)

    # namespaced => non-namespaced pairs
    pairs = [Pair(t...) for t in zip(all_spaced, all)]
    # filter out values which should be skipped (i.e. connected inputs)
    filter!(x -> first(x) ∉ skip, pairs)
    @show pairs

    # find duplicates and add namepsaces
    duplicates = Set()
    for i in 1:length(pairs)
        if last(pairs[i]) ∈ Set(last.(pairs[i+1:end]))
            push!(duplicates, last(pairs[i]))
        end
    end

    for i in 1:length(pairs)
        if pairs[i].second ∈ duplicates
            pairs[i] = pairs[i].first => pairs[i].first
        end
    end

    return pairs
end

end
