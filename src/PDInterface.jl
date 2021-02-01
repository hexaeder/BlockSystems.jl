module PDInterface

using IOSystems
using IOSystems: @check
using ModelingToolkit
using ModelingToolkit: getname, value
using NetworkDynamics
using PowerDynamics

export IONode, cconstruct_vertex

struct IONode <: AbstractNode
    iob::IOBlock
    dimensions::Int
    states::Vector
    params::Dict
end

function IONode(iob::IOBlock, params)
    states = iob.istates ∪ iob.outputs
    @parameters t i(t)
    @variables u_r(t) u_i(t)
    @check Set(i) == Set(iob.inputs) "inputs == i(t)"
    @check Set([u_r, u_i]) ⊆ Set(states) "IOBlock needs to contain output/istates u_r(t) u_i(t)"

    @check Set(keys(params)) == Set(iob.iparams) "params dict need to contain all iparams"

    states = vcat(value.([u_r, u_i]), states) |> unique
    inputs = vcat([value(i)], iob.inputs) |> unique

    IONode(iob, length(states), states, params)
end

function PowerDynamics.construct_vertex(ion::IONode)
    @parameters t i(t)
    # do we need a let block?
    # prop. only if we construct multiple vertices from the sam IONode
    first_states = ion.states
    first_inputs = [i]
    first_params = keys(ion.params) |> collect
    p_val = values(ion.params) |> collect

    gen = generate_io_function(ion.iob,
                               f_states=first_states,
                               f_inputs=first_inputs,
                               f_params=first_params)

    function rhs!(dx, x, e_s, e_d, p, t)
        i = total_current(e_s, e_d)
        gen.f_ip(dx, x, (i), p_val, t)
    end

    ODEVertex(f! = rhs!, dim=ion.dimensions, mass_matrix=gen.massm, sym=symbolsof(ion))
end

symbolsof(ion::IONode) = getname.(ion.states)
dimensions(ion::IONode) = ion.dimensions

end
