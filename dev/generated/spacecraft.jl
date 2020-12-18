using IOSystems
using ModelingToolkit
@parameters t, m
@variables x(t) F(t) v(t)
@derivatives D'~t

spacecraft = IOBlock([D(v) ~ F/m, D(x) ~ v], # define the equation
                     [F], # inputs of the system
                     [x], # outputs of the system
                     name = :spacecraft)
@show spacecraft.inputs
@show spacecraft.outputs
@show spacecraft.istates
nothing # hide

@parameters K
@variables o(t) i(t)

prop = IOBlock([o ~ K*i], [i], [o], name = :prop)

@variables Δ(t) p(t) m(t)
diff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)

propc = IOSystem([prop.i => diff.Δ], # conncect inptup of prop to output of diff
                 [diff, prop], # define subsystems
                 name=:propc)
propc.inputs_map
propc.istates_map
propc.outputs_map

a=[x,i,o]
ModelingToolkit.renamespace(:blubb,a[1])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

