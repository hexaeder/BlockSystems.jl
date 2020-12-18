#=
## Controlsystem: Spaceship
As an example we want to build an IOSystem controlling the altitude of a small spaccraft.
The spacecraft has mass ``m`` and can be controlled with thrusters wich applay the force ``F(t)`` to the spacecraft. The altitutude ``x(t)``
```math
\dot v(t) = \frac{ F(t) }{m}\\
\dot x(t) =  v(t)
```
in our model this system has the input ``F(t)``, the internal state ``v(t)`` and the output ``x(t)``.
=#

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

#=
We want to model a controler which takes a desired altituted as an input parameter and outputs the force for thrusters.

## Simple proportional controller
A proportional controller takes an input ``i`` and calculates the output proportional to the input.
```math
 o(t) = K\cdot i(t)
```
=#
@parameters K
@variables o(t) i(t)

prop = IOBlock([o ~ K*i], [i], [o], name = :prop)
#=
In order to make this usefull as an controller the input has to be the difference
between the reference and the system variable (negativ feedback). We can model this
as an IOSystem where
```math
Δ = p - m
```
=#
@variables Δ(t) p(t) m(t)
diff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)

#=
now we can connect both of the defined models to create an proportional controller
=#
propc = IOSystem([prop.i => diff.Δ], # conncect inptup of prop to output of diff
                 [diff, prop], # define subsystems
                 name=:propc)
propc.inputs_map
propc.istates_map
propc.outputs_map

a=[x,i,o]
ModelingToolkit.renamespace(:blubb,a[1])
