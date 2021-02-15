```@meta
CurrentModule = BlockSystems
```

# BlockSystems
## Basics
An input-output-system is characterized by a set of equations. These equations can be
either first order ODEs or explicit algebraic equations.
```math
\begin{aligned}
\dot \mathbf x(t) &= f(\mathbf x(t), \mathbf y(t), \mathbf i(t), p)\\
\mathbf y(t) &= g(\mathbf x(t), \mathbf y(t), \mathbf i(t), p)
\end{aligned}
```
such system contains of 
- states ($x$ and $y$) and
- parameters ($i$, $p$).

States are determined by the given equations (i.e. the equations describe how
the states change). Parameters are externally given. For IO systems we define subgroups
- states
  - internal states (`istates`) which are meant for internal use
  - output states (`outputs`) which might be used as inputs for other systems
- parameters
  - internal parameters (`iparams`) which are typically constant and
  - inputs (`inputs`) which can be connected to the outputs of other systems.
  
## Types
The base type for the BlockSystems is `AbstractIOSystem` with the has two concrete implementations. 
```@docs
AbstractIOSystem
```
An `IOBlock` consists of a set of equations, a set of inputs and outputs and a name.

```@docs
IOBlock
IOBlock(eqs::Vector{<:Equation}, inputs, outputs; name = gensym(:IOBlock))
```
An `IOSystem` consists of multiple `AbstractIOSystems` and the connections between them.
```@docs
IOSystem
IOSystem(cons, io_systems::Vector{<:AbstractIOSystem}; inputs_map = nothing, iparams_map = nothing, istates_map = nothing, outputs_map = nothing, name = gensym(:IOSystem), autopromote = true)
```

## Transformations
Each `IOSystem` can be transformed into an `IOBlock`. At this step the actual manipulation of the equations happen.
```@docs
connect_system
```

## Function building
```@docs
generate_io_function
```

## Block specifications
Sometimes it is useful to define a the input/output structure for a block.
```@docs
BlockSpec
fulfills
```

