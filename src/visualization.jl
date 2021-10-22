####
#### Pretty Printing
####
"""
Helper function for `Base.show` for IOBlock and IOSystem
"""
function print_variables(io::IO, V)
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
        diffs = rhs_differentials(iob)
        if !isempty(diffs)
            print(io, "\n  RHS differentials: "); print_variables(io, diffs)
        end
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
