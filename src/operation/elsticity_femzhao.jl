
using Revise, ApproxOperator, BenchmarkTools, Printf, SparseArrays
include("importmshwave.jl")
elements,nodes = import_fem("./msh/testzhao.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])
set𝝭!.(elements["Γᵗ"])
E = 2.1e6
ν=0.0
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)

prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->n₁*n₁)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->n₁*n₂)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->n₂*n₂)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->5.0)                 
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)    

ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
]
# m = zeros(2*nₚ,2*nₚ)
# kα = zeros(2*nₚ,2*nₚ)
# fα = zeros(2*nₚ)
#k = spzeros(2*nₚ,2*nₚ)
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)


ops[1](elements["Ω"],k)
ops[2](elements["Γ"],k,f)
ops[3](elements["Γᵗ"],f) 

d=k\f
d₁ = d[1:2:2*nₚ]
d₂ = d[2:2:2*nₚ] 
println("d₁ = ", d₁)
println("d₂ = ", d₂)
for ap in elements["Ω"]
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for (i,ξ) in enumerate(𝓖)
        if i == 1
            B₁ = ξ[:∂𝝭∂x]
            B₂ = ξ[:∂𝝭∂y]
            ε₁₁ = 0.0
            ε₂₂ = 0.0
            ε₁₂ = 0.0
            for (j,xⱼ) in enumerate(𝓒)
                ε₁₁ += B₁[j]*xⱼ.d₁
                ε₂₂ += B₂[j]*xⱼ.d₂
                ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
            end
            σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
            σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
            σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
            println("d₁ = ", d₁)
            println("d₂ = ", d₂)
            break
        end
    end
end
 


  