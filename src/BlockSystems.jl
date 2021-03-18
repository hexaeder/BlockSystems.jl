module BlockSystems

using LinearAlgebra
using DocStringExtensions
using ModelingToolkit
using ModelingToolkit: Parameter, ODESystem, Differential
using ModelingToolkit: rename, getname, renamespace, namespace_equation, namespace_equations, value, makesym, vars
using ModelingToolkit: equation_dependencies, asgraph, variable_dependencies, eqeq_dependencies, varvar_dependencies
using SymbolicUtils: Symbolic, operation
using LightGraphs

export AbstractIOSystem, IOBlock, IOSystem, independent_variable

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
        all_vars = vcat(inputs, outputs, istates, iparams, odes.iv)
        @check uniquenames(all_vars) "There seem to be a name clashes between inputs, iparams istates and outputs!"
        all_vars = Set(all_vars)
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
    map(eq->namespace_equation(eq, iob.name, independent_variable(iob)), eqs)
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

independent_variable(sys::IOSystem) = independent_variable(first(sys.systems))

"""
$(SIGNATURES)

Construct a new IOSystem from various subsystems.
Parameters
 - `cons`: the connections in the form `sub1.output => sub2.input`
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

    ivs = unique([independent_variable(sys) for sys in io_systems])
    length(ivs) == 1 || throw(ArgumentError("Multiple independent variables!"))

    nspcd_inputs = vcat([namespace_inputs(sys) for sys in io_systems]...)
    nspcd_iparams = vcat([namespace_iparams(sys) for sys in io_systems]...)
    nspcd_istates = vcat([namespace_istates(sys) for sys in io_systems]...)
    nspcd_outputs = vcat([namespace_outputs(sys) for sys in io_systems]...)
    nspcd_rem_states = vcat([namespace_rem_states(sys) for sys in io_systems]...)

    # check validity of provided connections
    @check Set(first.(cons)) ⊆ Set(nspcd_outputs) "First argument in connection needs to be output of subsystem."
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

Block specification, defines which inputs/outputs an [`AbstractIOSystem`](@ref)
should have. Contains two vectors of `Symbols`. Can be initialized with Vectors
of `Symbols`, `Num` or `<:Symbolic`.

Object is functor: call `(::BlockSpec)(ios)` to check wether `ios` fulfills
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
end
BlockSpec(in::Vector{Num}, out::Vector{Num}) = BlockSpec(value.(in), value.(out))
BlockSpec(in::Vector{<:Symbolic}, out::Vector{<:Symbolic}) = BlockSpec(getname.(in), getname.(out))

(bs::BlockSpec)(io) = fulfills(io, bs)

"""
    fulfills(io, bs::BlockSpec)::Bool

Check whether `io` fulfills the given [`BlockSpec`](@ref).
"""
fulfills(io, bs::BlockSpec) = Set(bs.inputs) ⊆ Set(getname.(io.inputs)) &&
                                Set(bs.outputs) ⊆ Set(getname.(io.outputs))


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


"""
Helper function for `Base.show` for IOBlock and IOSystem
"""
function print_variables(io::IO, V::Vector{Symbolic})
    printmax = get(io, :compact, false) ? 3 : 7
    l = length(V)
    if l == 0
        print(io, "(empty)")
        return
    elseif l > printmax
        print(io, V[1], ", …$(l-2)…,", V[end])
    else
        for (i, v) in enumerate(V)
            print(io, v)
            i < l && print(io, ", ")
        end
    end
end

function Base.show(io::IO, iob::IOBlock)
    compact = get(io, :compact, false)
    ioc = IOContext(io, :compact => true)
    if ~compact
        eqs = equations(iob.system)
        Base.printstyled(io, "IOBlock :$(iob.name) with $(length(eqs)) eqs", bold=true)
        print(io, "\n  ├ inputs:  "); print_variables(io, iob.inputs)
        print(io, "\n  ├ outputs: "); print_variables(io, iob.outputs)
        print(io, "\n  ├ istates: "); print_variables(io, iob.istates)
        print(io, "\n  └ iparams: "); print_variables(io, iob.iparams)
    else
        Base.printstyled(io, "$(iob.name): ", bold=true)
        print_variables(ioc, iob.inputs)
        print(io, " ↦ ")
        print_variables(ioc, iob.outputs)
    end
end

function Base.show(io::IO, ios::IOSystem)
    compact = get(io, :compact, false)
    ioc = IOContext(io, :compact => true)
    if ~compact
        Base.printstyled(io, "IOSystem :$(ios.name) with $(length(ios.systems)) subsystems", bold=true)
        for (i, sub) in enumerate(ios.systems)
            s = i == length(ios.systems) ? "└ " : "├ "
            print(ioc, "\n  ", s, sub)
        end
        Base.printstyled(io, "\nand $(length(ios.connections)) connections:", bold=true)
        for (i, con) in enumerate(ios.connections)
            s = i == length(ios.connections) ? "└ " : "├ "
            print(ioc, "\n  ", s, con.first, " ⇒ ", con.second)
        end
        all_unpromoted = Set(keys(filter(p->isequal(p.first, p.second), ios.namespace_map)))
        function show_promotions(vec)
            if isempty(vec)
                print(ioc, " (empty)")
                return
            end
            unpromoted = filter(v -> v ∈ all_unpromoted, vec)
            promoted = filter(v -> v ∉ all_unpromoted, vec)

            if ~isempty(unpromoted)
                s = isempty(promoted) ? "└ " : "├ "
                print(ioc, "\n  ", s)
                print_variables(io, unpromoted)
                print(ioc, " (unpromoted)")
            end

            for (i, v) in enumerate(promoted)
                s = i == length(promoted) ? "└ " : "├ "
                key = findfirst(val->isequal(val, v), ios.namespace_map)
                print(ioc, "\n  ", s, key, " ⇒ ", ios.namespace_map[key])
            end
        end
        Base.printstyled(io, "\ninputs:", bold=true)
        show_promotions(ios.inputs)
        Base.printstyled(io, "\noutputs:", bold=true)
        show_promotions(ios.outputs)
        Base.printstyled(io, "\nistates:", bold=true)
        show_promotions(ios.istates)
        Base.printstyled(io, "\niparams:", bold=true)
        show_promotions(ios.iparams)
    else
        Base.printstyled(io, "$(ios.name) (IOSystem): ", bold=true)
        print_variables(ioc, ios.inputs)
        print(io, " ↦ ")
        print_variables(ioc, ios.outputs)
    end
end

include("transformations.jl")
include("function_generation.jl")

end
