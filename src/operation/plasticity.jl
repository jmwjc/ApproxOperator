
"""
1D Plasticity
"""
function (op::Operator{:∫vₓσdx})(ap::T,k::AbstractMatrix{Float64},fint::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    K = op.K
    σy = op.σy
    tol = op.tol 
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        σₙ = ξ.σₙ
        αₙ = ξ.αₙ
        εᵖₙ = ξ.εᵖₙ
        𝑤 = ξ.𝑤
        Δεₙ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δεₙ += B[i]*xᵢ.Δd
        end
        # predict phase
        σᵗʳ = σₙ+E*Δεₙ
        fᵗʳ = abs(σᵗʳ) - (σy+K*αₙ)
        if fᵗʳ > tol
            Δγ = fᵗʳ/(E+K)
            ξ.σₙ = σᵗʳ - Δγ*E*sign(σᵗʳ)
            ξ.εᵖₙ = εᵖₙ + Δγ*sign(σᵗʳ)
            ξ.αₙ = αₙ + Δγ
            Eₜ = (E*K)/(E+K)
        else
            ξ.σₙ = σᵗʳ
            Eₜ = E
        end
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*Eₜ*B[j]*𝑤
            end
            fint[I] += B[i]*ξ.σₙ*𝑤
        end
    end
end

"""
morh-coulbom
"""
function (op::Operator{:∫vᵢσdΩ_mohr_coulomb})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    𝜙 = op.𝜙
    tol = op.tol 
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    C₁₂₁₂ =  μ
    Cᵗ₁₂₁₂ =  μ
    C₁₁₁₁ = Cᵢᵢᵢᵢ
    C₂₂₂₂ = Cᵢᵢᵢᵢ
    C₃₃₃₃ = Cᵢᵢᵢᵢ
    C₁₁₂₂ = Cᵢᵢⱼⱼ
    C₁₁₃₃ = Cᵢᵢⱼⱼ
    C₂₂₃₃ = Cᵢᵢⱼⱼ
    
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₃₃ = ξ.σ₃₃
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += 0.5*(B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁)
        end 
        # predict phase
        σ₁₁ᵗʳ =σ₁₁+C₁₁₁₁*Δε₁₁+C₁₁₂₂*Δε₂₂
        σ₂₂ᵗʳ =σ₂₂+C₂₂₂₂*Δε₂₂+C₁₁₂₂*Δε₁₁
        σ₁₂ᵗʳ =σ₁₂+C₁₂₁₂*Δε₁₂
        σ₃₃ᵗʳ =σ₃₃
        σ₁,σ₂,σ₃,n₁,n₂,n₃ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₃₃ᵗʳ,σ₁₂ᵗʳ,0.0,0.0)
        fᵗʳ = σ₁-σ₃ + (σ₁+σ₃) * sin(𝜙) - 2.0*c*cos(𝜙)
        if fᵗʳ > tol
            sin𝜙 = sin(𝜙)
            sin²𝜙 = sin(𝜙)*sin(𝜙)
            sin²𝜙psin𝜙 = sin²𝜙 +sin(𝜙) 
            sin²𝜙dsin𝜙 = sin²𝜙 -sin(𝜙) 
            ∂fC∂f=(4.0*sin²𝜙*λ+4.0*(1.0+sin²𝜙)*μ)
            Δγ = fᵗʳ/∂fC∂f
           
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*(2.0*(2.0*sin𝜙*λ+2.0*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)*n₃[1]*n₃[1])*μ))#上标1对应之前下标1
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*(2.0*(2.0*sin𝜙*λ+2.0*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)*n₃[2]*n₃[2])*μ))
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*(2.0*((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)*n₃[1]*n₃[2])*μ)
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*(((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)n₃[1]*n₃[1]))
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*(((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)n₃[2]*n₃[2]))
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*(((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)n₃[1]*n₃[2]))

            Cᵗ₁₁₁₁=C₁₁₁₁-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[1])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₃[1])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₃[1]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^4+4.0*(sin𝜙-1)^2*μ^2*(n₃[1])^4)/∂fC∂f
            Cᵗ₂₂₂₂=C₂₂₂₂-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[2])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₃[2])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₃[2]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^4+4.0*(sin𝜙-1)^2*μ^2*(n₃[2])^4)/∂fC∂f
            Cᵗ₁₁₂₂=C₁₁₂₂-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*((n₁[1])^2+(n₁[2])^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*((n₃[1])^2+(n₃[2])^2)+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^2*(n₁[2])^2-4.0*cos(𝜙)*cos(𝜙)*μ^2*((n₃[1])^2*(n₁[2])^2+(n₁[1])^2*(n₃[2])^2)+4.0*(sin𝜙-1)^2*μ^2*(n₃[1])^2*(n₃[2])^2)/∂fC∂f
        else
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = C₁₁₁₁
            Cᵗ₂₂₂₂ = C₂₂₂₂
            Cᵗ₁₁₂₂ = C₁₁₂₂
          
        end
        if isnan(σ₁₁)
           
              σ₁₁ = NaN
              error("程序终止：σ₁₁的值为NaN")
        end
          println(fint,k,σ₁₁)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₁₂₂*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₁₂₂*B₁[i]*B₁[j])*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
    
        end
    end
end    
