module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
import InteractiveUtils: subtypes

abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
# include("meshfree.jl")
# include("operation.jl")

include("preprocession/integration.jl")
include("preprocession/geometry.jl")
include("preprocession/importcomsol.jl")
include("approximation/tri3.jl")

export set𝝭!, set∇𝝭!

end