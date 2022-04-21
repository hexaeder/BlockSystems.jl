module BlockSystems

using LinearAlgebra
using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: ODESystem, Differential
using ModelingToolkit: get_iv, get_eqs, get_states
using ModelingToolkit: rename, getname, renamespace, namespace_equation, namespace_equations, value, vars
using ModelingToolkit: equation_dependencies, asgraph, variable_dependencies, eqeq_dependencies, varvar_dependencies
using ModelingToolkit.SymbolicUtils: Symbolic, operation, arguments, istree
using ModelingToolkit.Symbolics: tosymbol
using SciMLBase
using Graphs

export AbstractIOSystem, IOBlock, IOSystem, get_iv, equations

include("utils.jl")

"""
abstract supertype for [`IOBlock`](@ref) and [`IOSystem`](@ref).
"""
abstract type AbstractIOSystem end

function Base.getproperty(sys::AbstractIOSystem, name::Symbol)
    if name ∈ fieldnames(typeof(sys))
        return getfield(sys, name)
    end
    syms = vcat(getfield(sys, :inputs),
                getfield(sys, :iparams),
                getfield(sys, :istates),
                getfield(sys, :outputs),
                getfield(sys, :removed_states))

    i = findfirst(x -> getname(x) == name, syms)
    if i !== nothing
        return rename(syms[i], renamespace(getfield(sys, :name), name))
    end
    throw(error("Variable $name does not exist"))
end

function Base.getproperty(sys::AbstractIOSystem, sym::Symbolic)
    sym = remove_namespace(sys.name, sym)
    getproperty(sys, getname(sym))
end

Base.getproperty(sys::AbstractIOSystem, num::Num) = getproperty(sys, value(num))

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

    function IOBlock(name, inputs, iparams, istates, outputs, odes, rem_eqs; warn)
        @check name == getname(odes) "Name of inner ODESystem does not match name of IOBlock"
        @checkwarn warn Set(inputs) ⊆ Set(parameters(odes)) "Inputs should be parameters. You may ignore this warning if you want to specify an input which is not used in the eqs."
        @checkwarn warn Set(outputs) ⊆ Set(states(odes)) "Outputs should be variables. You may ignore this waring if you want to specify an output which is not used in the eqs (completly implicit)."
        @checkwarn warn Set(iparams) ⊆ Set(parameters(odes)) "iparams should be parameters of the eqs system"
        @checkwarn warn Set(istates) ⊆ Set(states(odes)) "istates should be variables of the eqs system"
        @check Set(parameters(odes)) ⊆ Set(inputs ∪ iparams) "params(eqs) ⊆ inputs ∪ iparams"
        @check Set(states(odes)) ⊆ Set(outputs ∪ istates) "states(eqs) ⊆ outputs ∪ istates"

        # lhs of the equation should be either differential of state/output
        for eq in equations(odes)
            (type, lhs_var) = eq_type(eq)
            if type === :explicit_diffeq
                @check lhs_var ∈ Set(istates ∪ outputs) "$eq : Lhs variable of diff eqs should be istate or output."
            elseif type === :explicit_algebraic
                @check lhs_var ∈ Set(istates ∪ outputs) "$eq : Lhs variable of explicit algebraic eqs should be istate or output."
            elseif type === :implicit_algebraic
                @checkwarn warn isnothing(lhs_var) "$eq : Please provide constraints in form `0 ~ f(...)`."
            end
        end

        additional_vars = get_variables([eq.rhs for eq in rem_eqs])
        all_vars = vcat(inputs, outputs, istates, iparams, get_iv(odes))
        @check uniquenames(all_vars) "There seem to be a name clashes between inputs, iparams istates and outputs!"
        all_vars = Set(all_vars)
        @checkwarn warn additional_vars ⊆ all_vars "Removed eqs should not contain new variables. This may be a bug since this feature is not fully implemented."

        if isempty(rem_eqs)
            rem_states = Symbolic[]
        else
            rem_states = lhs_var.(rem_eqs)
            @checkwarn warn isempty(rem_states ∩ all_vars) "removed states should not appear in in/out/is/ip equations"
            rem_rhs_vars = union([get_variables(eq.rhs) for eq in rem_eqs]...)
            if !(rem_rhs_vars ⊆ all_vars)
                remaining = setdiff(rem_rhs_vars, all_vars)
                @warn "RHS variables of removed eqs are not a subset of in/out/is/ip. Maybe $remaining got totally removed during simplification? This will lead to errors if you wan't to build functions for removed states."
            end
        end

        nometadata = check_metadata(equations(odes))
        if !isempty(nometadata)
            @warn "Transformation for $name dropped metadata in equations. Please report issue! Syms without data: $nometadata"
        end
        nometadata = check_metadata(rem_eqs)
        if !isempty(nometadata)
            @warn "Transformation for $name dropped metadata in removed equations. Please report issue! Syms without data: $nometadata"
        end

        # TODO: check IOBlock assumptions in inner constructor
        # - each state is represented by one lhs (static or diff) or implicit algebraic equation?
        # - lhs only first order or algebraic
        new(name, inputs, iparams, istates, outputs, odes, rem_states, rem_eqs)
    end
end

ModelingToolkit.get_iv(block::IOBlock) = get_iv(block.system)
ModelingToolkit.equations(block::IOBlock) = equations(block.system)

function namespace_rem_eqs(iob::IOBlock)
    eqs = iob.removed_eqs
    isempty(eqs) && return Equation[]

    map(eq->namespace_equation(eq, iob.system), eqs)
end

"""
$(SIGNATURES)

Construct a new IOBlock for the given arguments.

```@example
using BlockSystems, ModelingToolkit
@parameters t i(t)
@variables x(t) o(t)
D = Differential(t)

iob = IOBlock([D(x) ~ i, o ~ x], [i], [o], name=:iob)
```
"""
function IOBlock(eqs::Vector{<:Equation}, inputs, outputs; name = gensym(:IOBlock), iv = nothing, warn=true)
    IOBlock(name, eqs, inputs, outputs, Equation[]; iv, warn)
end

function IOBlock(name, eqs, inputs, outputs, rem_eqs; iv=nothing, warn=true)
    os = ODESystem(eqs, iv; name = name)

    inputs = value.(inputs) # gets the inputs as `tern` type
    outputs = value.(outputs) # gets the outputs as `tern` type
    istates = setdiff(get_states(os), outputs)
    iparams = setdiff(parameters(os), inputs)

    IOBlock(name, inputs, iparams, istates, outputs, os, rem_eqs; warn)
end


"""
$(SIGNATURES)

Construct a new IOBlock based on an existing. Deep-copy all fields and assigns new name.
"""
function IOBlock(iob::IOBlock; name=gensym(:IOBlock), warn=true)
    cp = deepcopy(iob)
    odes = ODESystem(get_eqs(cp.system), get_iv(cp), name=name)
    IOBlock(name, cp.inputs, cp.iparams, cp.istates, cp.outputs, odes, cp.removed_eqs; warn)
end

"""
$(TYPEDEF)

A composite `IOSystem` which consists of multiple [`AbstractIOSystem`](@ref) which are connected via
a vector of namespaced pairs (`subsys1.out => subsys2.in`).

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
    connections::Vector{Pair{Symbolic, Symbolic}}
    namespace_map::Dict{Symbolic, Symbolic}
    systems::Vector{AbstractIOSystem}
end
# TODO: check IOSystem assumptions in inner constructor?
# - check that every variable found in the equations is referenced by namespace map

ModelingToolkit.get_iv(sys::IOSystem) = get_iv(first(sys.systems))

"""
$(SIGNATURES)

Construct a new IOSystem from various subsystems.
Arguments:
 - `cons`:

    The connections in the form `sub1.output => sub2.input`. It is also possible to use simple
    algebraic equations such as `sub1.o1 + sub2.o2 => sub3.input`.

 - `io_systems`: Vector of subsystems
 - `namespace_map`: Provide collection of custom namespace promotions / renamings
   i.e. sub1.input => voltage`. Variables without entry in the map will be
   promoted automatically. Automatic promotion means that the sub-namespace is
   removed whenever it is possible without naming conflicts. The map may contain
   inputs, outputs, istates, iparams and removed_states. The rhs of the map can be provided
   as as Symbol: `sub1.input => :newname`.
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
                  autopromote = true)
    namespaces = [sys.name for sys in io_systems]
    allunique(namespaces) || throw(ArgumentError("Namespace collision in subsystems!"))

    ivs = unique([get_iv(sys) for sys in io_systems])
    length(ivs) == 1 || throw(ArgumentError("Multiple independent variables!"))

    nspcd_inputs = vcat([namespace_inputs(sys) for sys in io_systems]...)
    nspcd_iparams = vcat([namespace_iparams(sys) for sys in io_systems]...)
    nspcd_istates = vcat([namespace_istates(sys) for sys in io_systems]...)
    nspcd_outputs = vcat([namespace_outputs(sys) for sys in io_systems]...)
    nspcd_rem_states = vcat([namespace_rem_states(sys) for sys in io_systems]...)

    # check validity of provided connections
    @check vars(first.(cons)) ⊆ Set(nspcd_outputs) "First argument in connection may only contain outputs of subsystem."
    @check Set(last.(cons)) ⊆ Set(nspcd_inputs) "Second argument in connection needs to be input of subsystem."
    @check allunique(last.(cons)) "Multiple connections to same input!"

    # reduce the inputs to the open inputs
    open_inputs = setdiff(nspcd_inputs, last.(cons))

    # check validity of provided namespace maps
    namespace_map = fix_map_types(namespace_map)
    allowed_lhs = Set(open_inputs ∪ nspcd_iparams ∪ nspcd_istates ∪ nspcd_outputs ∪ nspcd_rem_states)
    @check keys(namespace_map) ⊆ Set(allowed_lhs) "namespace_map keys contain unknown symbols"
    @check uniquenames(values(namespace_map)) "naming conflict in rhs of user provided namespace promotions"

    # if the user provided a list of outputs, all other outputs become istates
    if outputs != :all
        outputs::Vector{Union{Symbol, Symbolic}} = value.(outputs)
        # check if outputs ∈ nspcd_outputs or referenced in rhs of namespace_map
        for (i, o) in enumerate(outputs)
            if o ∉ Set(nspcd_outputs)
                if o isa Symbol
                    key = findfirst(v->isequal(getname(v), o), namespace_map)
                else
                    key = findfirst(v->isequal(v, o), namespace_map)
                end
                # check if key references a output (namespace_map contains all)
                key ∉ Set(nspcd_outputs) && throw(ArgumentError("output $o ∉ outputs ∪ outputs_promoted"))
                outputs[i] = key
            end
        end
        former_outputs = setdiff(nspcd_outputs, outputs)
        nspcd_istates  = vcat(nspcd_istates, former_outputs)
        nspcd_outputs  = outputs |> collect
    end

    # auto promotion of the rest
    all_symbols = vcat(open_inputs, nspcd_iparams, nspcd_istates, nspcd_outputs, nspcd_rem_states)
    left_symbols = setdiff(all_symbols, keys(namespace_map))

    auto_promotions = create_namespace_promotions(collect(left_symbols), values(namespace_map), autopromote)

    namespace_map = merge(namespace_map, auto_promotions)

    # assert that there are no namespace clashes. this should be allways true!
    @assert uniquenames(keys(namespace_map)) "ERROR: complete namespace map lhs not unique"
    @assert uniquenames(values(namespace_map)) "ERROR: complete namespace map lhs not unique"

    nmspc_promote = ModelingToolkit.substituter(namespace_map)
    inputs  = nmspc_promote.(open_inputs)
    iparams = nmspc_promote.(nspcd_iparams)
    istates = nmspc_promote.(nspcd_istates)
    outputs = nmspc_promote.(nspcd_outputs)
    rem_states = nmspc_promote.(nspcd_rem_states)

    IOSystem(name,
             inputs, iparams, istates, outputs, rem_states,
             cons,
             namespace_map,
             io_systems)
end


export BlockSpec, fulfills

"""
    struct BlockSpec
    BlockSpec(in::Vector, out::Vector; in_strict=true, out_strict=false)

Block specification, defines which inputs/outputs an [`AbstractIOSystem`](@ref)
should have. Contains two vectors of `Symbols`. Can be initialized with Vectors
of `Symbols`, `Num` or `<:Symbolic`.

If `strict=true` the in/outputs *must equal* the specification. If `strict=false`
the block *must contain* the in/outputs from the specification.

Object is functor: call `(::BlockSpec)(ios)` to check whether `ios` fulfills
specification. See also [`fulfills`](@ref).

```
iob = IOBlock(...)
spec = BlockSpec([:uᵣ, :uᵢ], [:iᵣ, :iᵢ])
fulfills(iob, spec)
spec(iob)
```
"""
struct BlockSpec
    inputs::Vector{Symbol}
    outputs::Vector{Symbol}
    in_strict::Bool
    out_strict::Bool
end

BlockSpec(in::Vector, out::Vector; in_strict=true, out_strict=false) = BlockSpec(in, out, in_strict, out_strict)
BlockSpec(in::Vector{Num}, out::Vector{Num}; kwargs...) = BlockSpec(value.(in), value.(out); kwargs...)
BlockSpec(in::Vector{<:Symbolic}, out::Vector{<:Symbolic}; kwargs...) = BlockSpec(getname.(in), getname.(out); kwargs...)

(bs::BlockSpec)(io) = fulfills(io, bs)

"""
    fulfills(io, bs::BlockSpec)::Bool

Check whether `io` fulfills the given [`BlockSpec`](@ref).
"""
function fulfills(io, bs::BlockSpec)
    if bs.in_strict
        Set(bs.inputs) == Set(getname.(io.inputs))
    else
        Set(bs.inputs) ⊆ Set(getname.(io.inputs))
    end && if bs.out_strict
        Set(bs.outputs) == Set(getname.(io.outputs))
    else
        Set(bs.outputs) ⊆ Set(getname.(io.outputs))
    end
end


"""
    fix_map_types(map)

Creates Dict from `map`. Changes `Num`-types to `Symbolic`-types in the
user provided maps. Map can be Dict or  Array of Pairs.
If the rhs is of type `Symbol` it will be converted to a `Symbolic` type.
"""
function fix_map_types(map)
    dict = Dict(map)
    newdict = Dict{Symbolic, Symbolic}()
    for k in keys(dict)
        newkey = value(k)
        if dict[k] isa Symbol
            newval = rename(newkey, dict[k])
        else
            newval = value(dict[k])
        end
        newdict[newkey] = newval
    end
    newdict
end
fix_map_types(::Nothing) = Dict()

"""
    create_namepsace_promotions(syms, forbidden, autopromote)

Takes two lists of Symbols, `syms` and `forbidden`. Returns a Dict of namespace
promotions for each symbol in `syms`. Avoids name collisions inside the list as well
as with the `forbidden` symbols.
"""
function create_namespace_promotions(syms, forbidden, autopromote)
    dict = Dict()
    if autopromote
        promoted = remove_namespace.(syms)
        for i in 1:length(syms)
            forbidden_names = getname.(forbidden) ∪ getname.(promoted[[1:i-1;i+1:end]])
            dict[syms[i]] = getname(promoted[i])∈forbidden_names ? syms[i] : promoted[i]
        end
    else
        for i in 1:length(syms)
            dict[syms[i]] = syms[i]
        end
    end
    return dict
end

include("algebraic_elimination.jl")
include("transformations.jl")
include("function_generation.jl")
include("visualization.jl")

end
