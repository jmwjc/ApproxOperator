module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, similar
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type Approximator end
abstract type FiniteElement<:Approximator end
abstract type ReproducingKernel{T<:AbstractNode,𝒑,𝑠,𝜙} <: Approximator end
abstract type SpatialPartition end

include("node.jl")
include("meshfree.jl")
include("approximation.jl")
include("operation.jl")
include("import.jl")
export Node, SNode
export Poi1, Seg2, Tri3, Quad, PoiN, SegN
export importmsh
export RegularGrid
export Operator, prescribe!, similar, glue
export set𝓖!
export set𝝭!, set∇𝝭!, set∇̃𝝭!, set𝒏!

end
