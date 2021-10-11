#=
## Control System: Spaceship
As an example we want to build an IOSystem controlling the altitude of a spacecraft.
The spacecraft has mass ``m`` and can be controlled with thrusters which apply the force ``F(t)`` to the spacecraft. The altitude ``x(t)``
```math
\dot v(t) = \frac{ F(t) }{m}\\
\dot x(t) =  v(t)
```
in our model this system has the input ``F(t)``, the internal state ``v(t)`` (vertical velocity) and the output ``x(t)``.

```
       +------------+
F(t) --| spacecraft |-- x(t)
       | m, v(t)    |
       +------------+

```
=#

using BlockSystems
using ModelingToolkit
@parameters t M F(t)
@variables x(t) v(t)
D = Differential(t)

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v], # define the equation
                     [F], # inputs of the system
                     [x], # outputs of the system
                     name = :spacecraft)

#=
We want to model a controller which takes a desired altitude as an input parameter and outputs the force for thrusters.

## Simple proportional controller
A proportional controller takes an input `i` and calculates the output proportional to the input.
```math
 o(t) = K\cdot i(t)
```
```
       +--------+
i(t) --| prop K |-- o(t)
       +--------+
```
=#
@parameters K i(t)
@variables o(t)

prop = IOBlock([o ~ K*i], [i], [o], name = :prop)
nothing # hide
#=
In order to make this useful as an controller, the input has to be the difference
between the reference and the system variable (negative feedback). We can model this
as an IOSystem where
```math
Δ = p - m\,.
```
```
       +--------+
p(t) --|  diff  |-- Δ(t)
m(t) --|        |
       +--------+
```
=#
@parameters p(t) m(t)
@variables Δ(t)
diff = IOBlock([Δ ~ p - m], [p, m], [Δ], name=:diff)
nothing # hide
#=
Now we can connect both of the defined models to create an proportional controller

```
             +----------------------------------------+
             | propc                                  |
             |         +--------+   +--------+        |
  target(t)--|--p(t) --|  diff  |---| prop K |--o(t)--|--o(t)
feedback(t)--|--m(t) --|        |   +--------+        |
             |         +--------+                     |
             +----------------------------------------+
```

If we don't provide additional information the system will try to promote all of the enclosed
variables to the new systemwide namespace.
=#
prop_c = IOSystem([diff.Δ => prop.i], # connect output of diff to input of prop
                  [diff, prop], # subsystems
                  name=:propc)

#=
For finer control, it is often preferred to give new names manually, this is
done with the `namespace_map` argument. Per default, all of the outputs of the
subsystems will become outputs of the connected system (in this case also the
output `diff.Δ`). We can prevent this by supplying the `outputs` argument
manually. Sub outputs which are not referenced here will become internal states
of the connected system.

The rhs of the namespace map can be given as a Variable/Parameter type from MTK.
For simple renaming one can also give the rhs as a `Symbol` type.
=#
prop_c = IOSystem([diff.Δ => prop.i], [diff, prop],
                  namespace_map = [prop.o => o,
                                   diff.p => :target,
                                   diff.m => :feedback],
                  outputs = [o],
                  name=:propc)

#=
Right now, the created object is a container for the two included systems. However,
it is possible to transform the object into a new `IOBlock` by calling the `connect_system`
function. The resulting is equivalent to
```
             +----------------------------------+
  target(t)--| prop_c_block                     |--o(t)
feedback(t)--| o(t)=K*(target(t) - feedback(t)) |
             +----------------------------------+
```
=#
prop_c_block = connect_system(prop_c)
nothing #hide

#=
Now we can hook our spaceship to this controller. It does not matter whether we use the
connected `IOBlock` version `prop_c` or the `IOSystem` version `prop_c_block`. We want to build
the connected system

```
           +--------------------------------------------+
           |  control system                            |
           |       +--------+   +------------+          |
target(t)--|-------| prop_c |---| spacecraft |-x(t)--+--|--altitude(t)
           |  +-fb-|        |   | m, v(t)    |       |  |
           |  |    +--------+   +------------+       |  |
           |  +--------------------------------------+  |
           +--------------------------------------------+
```
=#
@variables altitude(t)
space_controller = IOSystem([prop_c.o => spacecraft.F, spacecraft.x => prop_c.feedback],
                            [prop_c, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])
## we want to reduce the space_controller to a block
space_controller = connect_system(space_controller)
@info "Variables of space_controller" space_controller equations(space_controller.system)

#=
## Simulate System
In order to simulate the system we can have to build the Julia functions.
=#
gen = generate_io_function(space_controller)
nothing # hide
#=
By doing so we get access to a named tuple with the fileds
- `gen.f_ip` in-place function
- `gen.f_oop` out-of-place function
- `gen.massm` mass matrix of the system
- `gen.states` symbols of states (in order)
- `gen.inputs` symbols of inputs (in order)
- `gen.params` symbols of parameters (in order)
- (see docstring for full list)

The functions have the form
`f_ip(du, u, inputs, params, t)`
where `u` are all the states (outputs stacked on top of internal states) and `t` is the independent variable of the system.
The order of the inputs and states can be controlled.
=#
gen = generate_io_function(space_controller, f_states=[altitude, v], f_params=[K, M])
@info "Generated function" gen.massm gen.states gen.inputs gen.params
nothing # hide

#=
Well, let's see how our model is doing.
=#
using Plots
using OrdinaryDiffEq
targetfun(t) = t>1.0 ? 1.0 : 0
odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [0.5, 1.0] # K, m
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(t->sol(t)[1],tspan..., label="altitude", title="proportional control")
plot!(t->targetfun(t),tspan..., label="target")

#=
Well who could have thought, proportional control looks like an harmonic oscillator 🤷‍♂️

## Defining a better controller
We might just add a damping term (a force proportional to the velocity of the spaceship).
If it works for a harmonic oscillator, it should work for our spaceship.
```
           +--------------------------------------------------------+
           |  control system                                        |
           |  +--------------------------------------------------+  |
           |  |    +--------+     +---+                          |  |
           |  +-v--| prop_v |-(-)-| d |     +------------+       |  |
           |       +--------+     | i |--F--| spacecraft |-v(t)--+  |
           |       +--------+     | f |     | m          |-x(t)--+--|--altitude(t)
target(t)--|-------| prop_c |-(+)-| f |     +------------+       |  |
           |  +-fb-|        |     +---+                          |  |
           |  |    +--------+                                    |  |
           |  +--------------------------------------------------+  |
           +--------------------------------------------------------+
```
In order to do so we have to slightly redefine the spaceship system: now the velocity `v(t)` is also an output and not and internal state.
=#

spacecraft = IOBlock([D(v) ~ F/M, D(x) ~ v],
                     [F],
                     [x,v],
                     name = :spacecraft)

# One can define new blocks based on previously defined blocks.
prop_v = IOBlock(prop, name=:prop_v)
fdiff = IOBlock(diff, name=:fdiff)

space_controller = IOSystem([spacecraft.v => prop_v.i,
                             spacecraft.x => prop_c.feedback,
                             prop_c.o => fdiff.p,
                             prop_v.o => fdiff.m,
                             fdiff.Δ => spacecraft.F],
                            [prop_v, prop_c, fdiff, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])

space_controller = connect_system(space_controller)
gen = generate_io_function(space_controller, f_states=[altitude, v], f_params=[prop_c.K, M, prop_v.K])

odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [1.0, 1.0, 0.5] # propc, M, propv
u0 = [0.0, 0.0] # altitude, v
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars=(0,1), label="altitude", title="better control")
plot!(t->targetfun(t),tspan..., label="target")

#=
## Defining an PT1 controller
```
             +-----------------------------------------------+
             | pi_c                                   +---+  |
             |                           +------------|   |  |
             |        +----+  +--------+ | +-------+  |sum|--|--o(t)
  target(t)--|--p(t)--|diff|--| prop K |-+-| int T |--|   |  |
feedback(t)--|--m(t)--|    |  +--------+   +-------+  +---+  |
             |        +----+                                 |
             +-----------------------------------------------+
```
=#

@parameters T a1(t) a2(t)
@variables Σ(t) altitude(t)
int = IOBlock([D(o) ~ 1/T * i - o], [i], [o], name=:int)
adder = IOBlock([Σ ~ a1 + a2], [a1, a2], [Σ], name=:add)

pi_c = IOSystem([diff.Δ => prop.i,
                 prop.o => int.i,
                 prop.o => adder.a1,
                 int.o => adder.a2],
                [diff, prop, int, adder],
                namespace_map = [diff.p => :target,
                                 diff.m => :feedback,
                                 adder.Σ => o],
                outputs = [o],
                name=:pi_c)
nothing # hide

# as before we can close the loop and build the control circuit

space_controller = IOSystem([pi_c.o => spacecraft.F , spacecraft.x => pi_c.feedback],
                            [pi_c, spacecraft],
                            namespace_map = [spacecraft.x => altitude],
                            outputs = [altitude])
space_controller = connect_system(space_controller, verbose=false)
@info "Variables of space_controller" space_controller equations(space_controller.system)

# and we can simulate and plot the system
gen = generate_io_function(space_controller, f_states=[altitude], f_params=[K, T, M])

odefun(du, u, p, t) = gen.f_ip(du, u, [targetfun(t)], p, t)
p = [0.5, -1.5, 1.0] # K, T, m
u0 = [0.0, 0.0, 0.0] # altitude, int.o, v
tspan = (0.0, 50.0)
prob = ODEProblem(odefun, u0, tspan, p)
sol = solve(prob, Tsit5())
plot(sol, vars=(0,[ 1,2 ]), label=["altitude" "integrator"], title="PT1 controller")
plot!(t->targetfun(t),tspan..., label="target")

# thank you for flying with us :)
