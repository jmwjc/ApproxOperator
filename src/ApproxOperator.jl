module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect, show
import InteractiveUtils: subtypes
import StaticArrays: SVector, SMatrix, @SArray
import Tensors: Vec, Tensor, SymmetricTensor, ⋅, ⊡, ⊗, gradient
import Printf: @printf
import Gmsh: gmsh

abstract type AbstractElement end
abstract type AbstractPiecewise<:AbstractElement end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
include("operation.jl")
include("meshfree.jl")

include("preprocession/integration.jl")
include("preprocession/geometry.jl")
include("preprocession/importcomsol.jl")
include("approximation/quad.jl")
include("approximation/quad8.jl")
include("approximation/tri3.jl")
include("approximation/tri6.jl")
include("approximation/tet4.jl")
include("approximation/hex8.jl")
include("approximation/seg2.jl")
include("approximation/seg3.jl")
include("approximation/poi1.jl")
include("approximation/reproducingkernel.jl")
include("approximation/kernelfunction.jl")
include("approximation/rkgradientsmoothing.jl")
include("approximation/Crouzeix-Raviart.jl")
include("approximation/grkgsi.jl")
include("approximation/rkwave.jl")
include("approximation/piecewise.jl")
include("operation/potential.jl")
include("operation/thin_shell.jl")
include("operation/elasticity.jl")
include("operation/incompressible.jl")
include("operation/curved beam.jl")
include("operation/plasticity.jl")
include("operation/hyperelasticity.jl")
include("operation/phasefield.jl")
include("operation/error_estimates.jl")
include("preprocession/importmsh.jl")
include("operation/Kirchhoff_Love_plate.jl")
include("operation/Mindlin_Reissner_plate.jl")
include("operation/plasticity.jl")
include("operation/thin_shell_penalty.jl")

export prescribe!
export Operator
export 𝑿ᵢ, 𝑿ₛ
export Element
export TRElement
export ReproducingKernel, RKGradientSmoothing, GRKGradientSmoothing, PiecewiseParametric, PiecewisePolynomial, RegularGrid
export set𝝭!, set∇𝝭!, set∇²𝝭!, set∇̂³𝝭!
export getPhysicalGroups, get𝑿ᵢ, getElements, addEdgeElements
export getPiecewiseElements, getPiecewiseBoundaryElements
export getMacroElements, getMacroBoundaryElements, getCurvedElements, getCurvedPiecewiseElements

end