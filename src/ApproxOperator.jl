module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect, show
import InteractiveUtils: subtypes
import Printf: @printf

abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
include("operation.jl")
include("meshfree.jl")

include("preprocession/integration.jl")
include("preprocession/geometry.jl")
include("preprocession/field.jl")
include("preprocession/importcomsol.jl")
include("approximation/quad.jl")
include("approximation/quad8.jl")
include("approximation/tri3.jl")
include("approximation/seg2.jl")
include("approximation/seg3.jl")
include("approximation/poi1.jl")
include("approximation/reproducingkernel.jl")
include("approximation/kernelfunction.jl")
include("approximation/rkgradientsmoothing.jl")
include("approximation/Crouzeix-Raviart.jl")
include("approximation/grkgsi.jl")
include("approximation/rkwave.jl")
include("operation/potential.jl")
include("operation/elasticity.jl")
include("operation/curved beam.jl")
include("operation/hyperelasticity.jl")
include("operation/error_estimates.jl")
include("preprocession/importmsh.jl")
include("operation/Kirchhoff_Love_plate.jl")

include("littletools.jl")


export prescribe!
export Operator
export Node, Element
export TRElement
export ReproducingKernel, RKGradientSmoothing, GRKGradientSmoothing
export set𝝭!, set∇𝝭!, set∇²𝝭!
export 𝑿ᵢ,𝑿ₛ
end