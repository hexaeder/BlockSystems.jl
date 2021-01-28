module IOSystems

using LinearAlgebra
using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: Parameter, ODESystem, Differential
using ModelingToolkit: rename, getname, renamespace, namespace_equation, namespace_equations, value, makesym, vars
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
    syms = vcat(sys.inputs, sys.iparams, sys.istates, sys.outputs, sys.removed_states)

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
namespace_rem_states(ios::AbstractIOSystem) = namespace_property(ios, :removed_states)

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
    removed_states::Vector{Symbolic}
    removed_eqs::Vector{Equation}

    function IOBlock(name, inputs, iparams, istates, outputs, odes, rem_eqs)
        @check name == odes.name "Name of inner ODESystem does not match name of IOBlock"
        @check Set(inputs) ⊆ Set(parameters(odes)) "inputs must be parameters"
        @check Set(outputs) ⊆ Set(states(odes)) "outputs must be variables"
        @check Set(iparams) ⊆ Set(parameters(odes)) "iparams must be parameters"
        @check Set(istates) ⊆ Set(states(odes)) "istates must be variables"
        @check Set(inputs ∪ iparams) == Set(parameters(odes)) "inputs ∪ iparams != params"
        @check Set(outputs ∪ istates) == Set(states(odes)) "outputs ∪ istates != states"

        additional_vars = vars([eq.rhs for eq in rem_eqs])
        all_vars = Set(inputs ∪ outputs ∪ istates ∪ iparams ∪ (odes.iv, ))
        @check additional_vars ⊆ all_vars "removed eqs should not contain new variables"

        if isempty(rem_eqs)
            rem_states = Symbolic[]
        else
            rem_states = lhs_var.(rem_eqs)
            @check isempty(rem_states ∩ all_vars) "removed states should not appear in in/out/is/ip"
            rem_rhs_vars = union([vars(eq.rhs) for eq in rem_eqs]...)
            @check rem_rhs_vars ⊆ all_vars "rhs of removed eqs should be subset of in/out/is/ip"
        end

        # TODO: check IOBlock assumptions in inner constructor
        # - each state is represented by one lhs (static or diff) or implicit algebraic equation?
        # - lhs only first order or algebraic
        new(name, inputs, iparams, istates, outputs, odes, rem_states, rem_eqs)
    end
end

independent_variable(block::IOBlock) = block.system.iv

function namespace_rem_eqs(iob::IOBlock)
    eqs = iob.removed_eqs
    isempty(eqs) && return Equation[]
    map(eq->namespace_equation(eq, iob.name, independent_variable(iob).name), eqs)
end

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
function IOBlock(eqs::Vector{<:Equation}, inputs, outputs; name = gensym(:IOBlock))
    IOBlock(name, eqs, inputs, outputs, Equation[])
end

function IOBlock(name, eqs, inputs, outputs, rem_eqs)
    os = ODESystem(eqs, name = name)

    inputs = value.(inputs) # gets the inputs as `tern` type
    outputs = value.(outputs) # gets the outputs as `tern` type
    istates = setdiff(os.states, outputs)
    iparams = setdiff(parameters(os), inputs)

    IOBlock(name, inputs, iparams, istates, outputs, os, rem_eqs)
end


"""
$(SIGNATURES)

Construct a new IOBlock based on an existing. Deep-copy all fields and assigns new name.
"""
function IOBlock(iob::IOBlock; name=gensym(:IOBlock))
    cp = deepcopy(iob)
    odes = ODESystem(cp.system.eqs, cp.system.iv, name=name)
    IOBlock(name, cp.inputs, cp.iparams, cp.istates, cp.outputs, odes, cp.removed_eqs)
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
    removed_states::Vector{Symbolic}
    connections::Dict{Symbolic, Symbolic}
    namespace_map::Dict{Symbolic, Symbolic}
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
 - `namespace_map`: Provide collection of custom namespace promotions / renamings
   i.e. `sub1.input => voltage`. Variables without entry in the map will be
   promoted automatically. Automatic promotion means that the sub-namespace is
   removed whenever it is possible without naming conflicts. The map may contain
   inputs, outputs, istates, iparams and removed_states. TODO: Allow :symbols in rhs.
 - `outputs`: Per default, all of the subsystem outputs will become system outputs. However,
   by providing a list of variables as outputs *only these* will become outputs of the
   new system. All other sub-outputs will become internal states of the connected system
   (and might be optimized away in `connect_system`).
 - `name`: namespace
 - `autopromote=true`: enable/disable automatic promotion of variable names to system namespace
"""
function IOSystem(cons,
                  io_systems::Vector{<:AbstractIOSystem};
                  namespace_map = nothing,
                  outputs = :all,
                  name = gensym(:IOSystem),
                  autopromote = true
                  )
    namespaces = [sys.name for sys in io_systems]
    @check allunique(namespaces) "Namespace collision in subsystems!"

    ivs = unique([independent_variable(sys) for sys in io_systems])
    @check length(ivs) == 1 "Multiple independent variables!"

    @check allunique(first.(cons)) "Multiple connections to same input!"
    nspcd_inputs = vcat([namespace_inputs(sys) for sys in io_systems]...)
    nspcd_iparams = vcat([namespace_iparams(sys) for sys in io_systems]...)
    nspcd_istates = vcat([namespace_istates(sys) for sys in io_systems]...)
    nspcd_outputs = vcat([namespace_outputs(sys) for sys in io_systems]...)
    nspcd_rem_states = vcat([namespace_rem_states(sys) for sys in io_systems]...)

    # check validity of provided connections
    @check Set(first.(cons)) ⊆ Set(nspcd_inputs) "First argument in connection needs to be input of subsystem."
    @check Set(last.(cons)) ⊆ Set(nspcd_outputs) "Second argument in connection needs to be output of subsystem."
    # reduce the inputs to the open inputs
    open_inputs = setdiff(nspcd_inputs, first.(cons))

    # check validity of provided namespace maps
    namespace_map = fix_map_types(namespace_map)
    allowed_lhs = Set(open_inputs ∪ nspcd_iparams ∪ nspcd_istates ∪ nspcd_outputs ∪ nspcd_rem_states)
    @check keys(namespace_map) ⊆ Set(allowed_lhs) "namespace_map keys contain unknown symbols"
    @check uniquenames(values(namespace_map)) "naming conflict in rhs of user provided namespace promotions"

    # if the user provided a list of outputs, all other outputs become istates
    if outputs != :all
        # check if outputs ∈ nspcd_outputs or referenced in rhs of namespace_map
        for (i, o) in enumerate(outputs)
            if o ∉ Set(nspcd_outputs)
                key = findfirst(v->isequal(v, o), namespace_map)
                @check key ∈ Set(nspcd_outputs) "Output $o ∉ outputs ∪ outputs_promoted"
                outputs[i] = key
            end
        end
        former_outputs = setdiff(nspcd_outputs, outputs)
        nspcd_istates  = vcat(nspcd_istates, former_outputs)
        nspcd_outputs  = outputs |> collect
    end

    # TODO continue here

    # auto promotion of the rest
    all_symbols = vcat(open_inputs, nspcd_iparams, nspcd_istates, nspcd_outputs, nspcd_rem_states)
    left_symbols = setdiff(all_symbols, keys(user_promotions))
    auto_promotions = create_namespace_promotions(collect(left_symbols), values(user_promotions))
    promotions = merge(user_promotions, auto_promotions)


    in_map = Dict()
    for s in open_inputs
        in_map[s] = promotions[s]
    end
    ip_map = Dict()
    for s in nspcd_iparams
        ip_map[s] = promotions[s]
    end
    is_map = Dict()
    for s in nspcd_istates
        is_map[s] = promotions[s]
    end
    out_map = Dict()
    for s in nspcd_outputs
        out_map[s] = promotions[s]
    end
    rem_map = Dict()
    for s in nspcd_rem_states
        rem_map[s] = promotions[s]
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
        collect(values(rem_map)),
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
