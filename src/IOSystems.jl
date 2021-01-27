module IOSystems

using LinearAlgebra
using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: Parameter, ODESystem, Differential
using ModelingToolkit: rename, getname, renamespace, namespace_equations, value, makesym, vars
using ModelingToolkit: equation_dependencies, asgraph, variable_dependencies, eqeq_dependencies, varvar_dependencies
using SymbolicUtils: Symbolic
using LightGraphs

export AbstractIOSystem, IOBlock, IOSystem

include("utils.jl")

"""
abstract supertype for [`IOBlock`](@ref) and [`IOSystem`](@ref).
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

A basic IOSystem which consists of a single ODESystem.
$(TYPEDFIELDS)
"""
struct IOBlock <: AbstractIOSystem
    name::Symbol
    inputs::Vector{Symbolic}
    iparams::Vector{Symbolic}
    istates::Vector{Symbolic}
    outputs::Vector{Symbolic}
    system::ODESystem
    removed_eqs::Vector{Equation}

    function IOBlock(name, inputs, iparams, istates, outputs, os, removed_eqs)
        @check Set(inputs) ⊆ Set(parameters(os)) "inputs must be parameters"
        @check Set(outputs) ⊆ Set(states(os)) "outputs must be variables"
        @check Set(iparams) ⊆ Set(parameters(os)) "iparams must be parameters"
        @check Set(istates) ⊆ Set(states(os)) "istates must be variables"
        @check Set(inputs ∪ iparams) == Set(parameters(os)) "inputs ∪ iparams != params"
        @check Set(outputs ∪ istates) == Set(states(os)) "outputs ∪ istates != states"

        additional_vars = vars([eq.rhs for eq in removed_eqs])
        all_vars = Set(inputs ∪ outputs ∪ istates ∪ iparams ∪ (os.iv, ))
        @check additional_vars ⊆ all_vars "removed eqs should not contain new variables"

        # TODO: check IOBlock assumptions in inner constructor
        # - each state is represented by one lhs (static or diff) or implicit algebraic equation?
        # - lhs only first order or algebraic
        new(name, inputs, iparams, istates, outputs, os, removed_eqs)
    end
end

independent_variable(block::IOBlock) = block.system.iv

"""
$(SIGNATURES)

Construct a new IOBlock for the given arguments.

```@example
using IOSystems, ModelingToolkit
@parameters t i(t)
@variables x(t) o(t)
@derivatives D'~t

iob = IOBlock([D(x) ~ i, o ~ x], [i], [o], name=:iob)
```
"""
function IOBlock(eqs::AbstractVector{<:Equation}, inputs, outputs;
                 name = gensym(:IOBlock), removed_eq = Equation[])
    os = ODESystem(eqs, name = name)

    inputs = value.(inputs) # gets the inputs as `tern` type
    outputs = value.(outputs) # gets the outputs as `tern` type
    istates = setdiff(os.states, outputs)
    iparams = setdiff(parameters(os), inputs)

    IOBlock(name, inputs, iparams, istates, outputs, os, removed_eq)
end

"""
$(SIGNATURES)

Construct a new IOBlock based on an existing. Deep-copy all fields and assigns new name.
"""
function IOBlock(iob::IOBlock; name=gensym(:IOBlock))
    cp = deepcopy(iob)
    IOBlock(cp.system.eqs, cp.inputs, cp.outputs, name=name)
end

"""
$(TYPEDEF)

A composite `IOSystem` which consists of multiple [`AbstractIOSystem`](@ref) which are connected via
a vector of namespaced pairs (`subsys1.out1 => subsys2.in1`).

An `IOSystem` contains maps how to promote the namespaced variables of the subsystem to the new scope
   subsys1₊x(t) => x(t)
   subsys1₊y(t) => subsys1₊y(t)
   subsys2₊y(t) => subsys2₊y(t)

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
# TODO: check IOSystem assumptions in inner constructor?
# - check that every variable found in the equations is referenced by namespace map

independent_variable(sys::IOSystem) = independent_variable(first(sys.systems))

"""
$(SIGNATURES)

Construct a new IOSystem from various subsystems.
Parameters
 - `cons`: the connections in the form `sub2.input => sub2.output`
 - `io_systems`: Vector of subsystems
 - `inputs_map`, `iparams_map`, `istates_map`:
   Provide custom namespace promotion / renaming i.e. `sub1.input => voltage`.
   Variables without entry in the map will be promoted automatically.
 - `outputs_map`:
   Provide custom namespace promotion for outputs. If this parameter is not given,
   all of the outputs (connected and unconnected) will be marked as outputs of the whole
   system. If the outputs map is given only the given outputs will be marked as outputs.
   All other outputs of subsystems will become internal states of the connected system
   (and might be optimized away in `connect_system`).
 - `name`: namespace
 - `autopromote=true`: enable/disable automatic promotion of variable names to system namespace
"""
function IOSystem(cons,
                  io_systems::Vector{<:AbstractIOSystem};
                  inputs_map = nothing,
                  iparams_map = nothing,
                  istates_map = nothing,
                  outputs_map = nothing,
                  name = gensym(:IOSystem),
                  autopromote = true
                  )
    namespaces = [sys.name for sys in io_systems]
    @check allunique(namespaces) "Namespace collision in subsystems!"

    ivs = unique([independent_variable(sys) for sys in io_systems])
    @check length(ivs) == 1 "Multiple independent variables!"

    @check allunique(first.(cons)) "Multiple connections to same input!"
    namespaced_inputs = vcat([namespace_inputs(sys) for sys in io_systems]...)
    namespaced_iparams = vcat([namespace_iparams(sys) for sys in io_systems]...)
    namespaced_istates = vcat([namespace_istates(sys) for sys in io_systems]...)
    namespaced_outputs = vcat([namespace_outputs(sys) for sys in io_systems]...)

    # check validity of provided connections
    @check Set(first.(cons)) ⊆ Set(namespaced_inputs) "First argument in connection needs to be input of subsystem."
    @check Set(last.(cons)) ⊆ Set(namespaced_outputs) "Second argument in connection needs to be output of subsystem."
    # reduce the inputs to the open inputs
    open_inputs = setdiff(namespaced_inputs, first.(cons))

    # check validity of provided namespace maps
    inputs_map = fix_map_types(inputs_map)
    @check keys(inputs_map) ⊆ Set(open_inputs) "inputs_map !⊆ open_inputs"
    iparams_map = fix_map_types(iparams_map)
    @check keys(iparams_map) ⊆ Set(namespaced_iparams) "iparams_map !⊆ iparams"
    istates_map = fix_map_types(istates_map)
    @check keys(istates_map) ⊆ Set(namespaced_istates) "istates_map !⊆ istates"
    outputs_map = fix_map_types(outputs_map)
    @check keys(outputs_map) ⊆ Set(namespaced_outputs) "outputs_map !⊆ outputs"
    user_promotions = merge(inputs_map, iparams_map, istates_map, outputs_map)
    @check uniquenames(values(user_promotions)) "naming conflict in rhs of user provided namespace promotions"

    # if the user provided a map of outputs, those outputs which are note referenced become states!
    if !isempty(outputs_map)
        former_outputs = setdiff(namespaced_outputs, keys(outputs_map))
        namespaced_istates = vcat(namespaced_istates, former_outputs)
        namespaced_outputs = keys(outputs_map) |> collect
    end

    # auto promotion of the rest
    if autopromote
        all_symbols = vcat(open_inputs, namespaced_iparams, namespaced_istates, namespaced_outputs)
        left_symbols = setdiff(all_symbols, keys(user_promotions))
        auto_promotions = create_namespace_promotions(collect(left_symbols), values(user_promotions))
        promotions = merge(user_promotions, auto_promotions)
    else
        promotions = user_promotions
    end

    in_map = Dict()
    for s in open_inputs
        in_map[s] = promotions[s]
    end
    ip_map = Dict()
    for s in namespaced_iparams
        ip_map[s] = promotions[s]
    end
    is_map = Dict()
    for s in namespaced_istates
        is_map[s] = promotions[s]
    end
    out_map = Dict()
    for s in namespaced_outputs
        out_map[s] = promotions[s]
    end

    # assert that there are no namespace clashes. this should be allways true!
    @assert uniquenames(values(inputs_map)) "namespace promotion of inputs clashed with manually given inputs_map"
    @assert uniquenames(values(iparams_map)) "namespace promotion of iparams clashed with manually given iparams_map"
    @assert uniquenames(values(istates_map)) "namespace promotion of istates clashed with manually given istates_map"
    @assert uniquenames(vcat(collect.(keys.([inputs_map, iparams_map, istates_map, outputs_map]))...)) "lhs of namespace promotion not unique"
    @assert uniquenames(vcat(collect.(values.([in_map, ip_map, is_map, out_map]))...)) "rhs of namespace promotion not unique"

    IOSystem(
        name,
        collect(values(in_map)),
        collect(values(ip_map)),
        collect(values(is_map)),
        collect(values(out_map)),
        Dict(cons),
        in_map,
        ip_map,
        is_map,
        out_map,
        io_systems,
      )
end

"""
    fix_map_types(map)

Creates Dict from `map`. Changes `Num`-types to `Symbolic`-types in the
user provided maps. Map can be Dict or  Array of Pairs.
"""
function fix_map_types(map)
    dict = Dict(map)
    newdict = Dict{Symbolic, Symbolic}()
    for k in keys(dict)
        newkey = value(k)
        newdict[newkey] = value(dict[k])
    end
    newdict
end
fix_map_types(::Nothing) = Dict()

"""
    create_namepsace_promotions(syms, forbidden)

Takes two lists of Symbols, `syms` and `forbidden`. Returns a Dict of namespace
promotions for each symbol in `syms`. Avoids name collisions inside the list as well
as with the `forbidden` symbols.
"""
function create_namespace_promotions(syms, forbidden)
    promoted = remove_namespace.(syms)
    dict = Dict()
    for i in 1:length(syms)
        forbidden_names = getname.(forbidden) ∪ getname.(promoted[[1:i-1;i+1:end]])
        dict[syms[i]] = getname(promoted[i])∈forbidden_names ? syms[i] : promoted[i]
    end
    return dict
end

include("transformations.jl")
include("function_generation.jl")

end
