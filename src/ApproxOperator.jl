module ApproxOperator

using StaticArrays
using SparseArrays
using LinearAlgebra

abstract type AbstractNode end
abstract type PhysicalNode<:AbstractNode end
abstract type ParametricNode<:AbstractNode end
abstract type Approximator end
abstract type BasisFunction end
abstract type KernelFunction end
abstract type ShapeFunction end
abstract type SpatialPartition end
abstract type AbstractFunction end
abstract type AbstractPool end
abstract type AbstractPoi <: Approximator end
abstract type AbstractSeg <: Approximator end
abstract type AbstractTri <: Approximator end
abstract type AbstractQuad <: Approximator end
abstract type Operator end

include("nodes.jl")
export set_nodes!
include("approximation.jl")
include("operation.jl")
export EBCDOFS, Potential_Ω, Potential_Γᵍ, Potential_Γᵗ, Potential_Γᵍ_Nitsche, Potential_Γᵍ_penalty, Potential_Γᵍ_Lagrange_multiplier, PlaneStress_Ω, PlaneStrain_Ωᵛ, PlaneStrain_Ωᵈ, PlaneStress_Γᵗ,PlaneStress_Γᵍ_penalty, NonlinearElasticity_Γᵍ_penalty, NonlinearPlaneStress_Γᵍ_penalty, NonlinearPlaneStress_C_Ω
export L₂Error_scale, H₁Error_scale, L₂Error_tensor, L₂Error_2nd_order_tensor, HₑError_PlaneStress
export PlaneStress_PhaseField_Ω, SecondOrderPhaseField, Update_HistoryField_PlaneStress, Update_Friction_PhaseField_PlaneStress
include("import.jl")
export import_msh
include("export.jl")
export export_shape_functions, VTKExport

end
