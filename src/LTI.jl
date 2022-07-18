export identify_lti

function _state_matrix(expr, vars)
    if eltype(expr) <: Equation
        expr = getproperty.(expr, :rhs)
    end

    A = Matrix(undef, length(expr), length(vars))
    for (ieq, ex) in enumerate(expr)
        for (istate, state) in enumerate(vars)
            coeff = (ex - substitute(ex, state => 0))/state |> simplify
            if !isempty(Set(get_variables(coeff)) ∩ Set(vars))
                error("Could not properly extract $state from $(ex), got $coeff")
            end
            A[ieq, istate] = coeff
        end
    end
    return narrow_type(A)
end


"""
    identify_lti(blk::IOBlock)

Identify the matrices A,B,C and D from IOBlock.

"""
function identify_lti(blk::IOBlock)
    inputs  = blk.inputs
    outputs = blk.outputs
    output_rhs = similar(outputs, Any)
    odestates = Term[]
    odeidx    = Int[]
    algstates = Term[]
    algidx    = Int[]
    for (i, eq) in enumerate(equations(blk))
        (type, var) = BlockSystems.eq_type(eq)
        if type == :explicit_diffeq
            push!(odestates, var)
            push!(odeidx, i)
        elseif type == :explicit_algebraic
            @assert var ∈ Set(outputs) "Explicit algebraic equations musst represent outputs!"
            push!(algstates, var)
            push!(algidx, i)
        else
            error("Only pure LTI systems. Can not handle $type-type equations!")
        end
        #check if variable is an output
        outidx = findfirst(isequal(var), outputs)
        if !isnothing(outidx) # is output variable
            if type == :explicit_diffeq
                output_rhs[outidx] = var
            elseif type == :explicit_algebraic
                output_rhs[outidx] = eq.rhs
            else
                error()
            end
        end
    end

    for eq in equations(blk)
        @assert isempty(Set(get_variables(eq.rhs)) ∩ algstates) "The rhs of the equations must not contain algebraic states!"
    end

    A = _state_matrix(equations(blk)[odeidx], odestates)
    B = _state_matrix(equations(blk)[odeidx], inputs)
    C = _state_matrix(output_rhs, odestates)
    D = _state_matrix(output_rhs, inputs)

    Avars = Set(mapreduce(get_variables, vcat, A))
    Bvars = Set(mapreduce(get_variables, vcat, B))
    Cvars = Set(mapreduce(get_variables, vcat, C))
    Dvars = Set(mapreduce(get_variables, vcat, D))

    @assert Avars ⊆ Set(blk.iparams) "Matrix A contains non-parameters. Thats an error!"
    @assert Bvars ⊆ Set(blk.iparams) "Matrix B contains non-parameters. Thats an error!"
    @assert Cvars ⊆ Set(blk.iparams) "Matrix C contains non-parameters. Thats an error!"
    @assert Dvars ⊆ Set(blk.iparams) "Matrix D contains non-parameters. Thats an error!"

    return A, B, C, D
end
