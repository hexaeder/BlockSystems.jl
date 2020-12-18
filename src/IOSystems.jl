module IOSystems

using DocStringExtensions
using ModelingToolkit
# import internal stuff from ModelingToolkit
import ModelingToolkit.rename, ModelingToolkit.getname, ModelingToolkit.renamespace

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
    syms = vcat(sys.inputs, sys.outputs, sys.istates)

    i = findfirst(x -> getname(x) == name, syms)
    if i !== nothing
        return rename(syms[i], renamespace(sys.name, name))
    end
    throw(error("Variable $name does not exist"))
end

namespace_inputs(ios::AbstractIOSystem) = [renamespace(ios.name,x) for x in ios.inputs]
namespace_istates(ios::AbstractIOSystem) = [renamespace(ios.name,x) for x in ios.istates]
namespace_outputs(ios::AbstractIOSystem) = [renamespace(ios.name,x) for x in ios.outputs]

"""
$(TYPEDEF)

A basic IOSystem which consists of a single ODESystems.
$(TYPEDFIELDS)
"""
struct IOBlock <: AbstractIOSystem
    name::Symbol
    inputs::Vector{Term}
    istates::Vector{Term}
    outputs::Vector{Term}
    system::ODESystem
end

function IOBlock(eqs::AbstractVector{<:Equation}, inputs, outputs; name = gensym(:IOBlock))
    os = ODESystem(eqs, name = name)

    # put asserts here

    inputs = os.states ∩ inputs # gets the inputs as `tern` type
    outputs = os.states ∩ outputs # gets the outputs as `tern` type
    istates = setdiff(os.states, inputs ∪ outputs)

    @assert Set(os.states) == Set(inputs ∪ outputs ∪ istates)

    IOBlock(name, inputs, istates, outputs, os)
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
    inputs::Vector{Term}
    istates::Vector{Term}
    outputs::Vector{Term}
    connections::Dict{Term, Term}
    inputs_map::Dict{Term, Term}
    istates_map::Dict{Term, Term}
    outputs_map::Dict{Term, Term}
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

function collect_and_resolve_namespace(io_systems, property; skip = [])
    all = vcat((getproperty(ios, property) for ios in io_systems)...)
    all_spaced = vcat((
        map(x -> renamespace(ios.name, x), getproperty(ios, property)) for
        ios in io_systems
    )...)


    l = length(all)
    filter!(x->x ∉ Set(skip), all)
    if skip != []
        @info skip all
    end
    if l - length(all) ≠ 0
        @info "gefiltert!"
    end

    duplicates = Set()
    for i = 1:length(all)
        if all[i] ∈ Set(all[i+1:end])
            push!(duplicates, all[i])
        end
    end
    map_resolved = Dict()
    for ios in io_systems
        vec = getproperty(ios, property)
        resolved = map(x -> x ∉ duplicates ? x : renamespace(ios.name, x), vec)
        namespaced = map(x -> renamespace(ios.name, x), vec)
        for kv in zip(resolved, namespaced)
            push!(map_resolved, Pair(kv...))
        end
    end
    map_resolved
end

end
