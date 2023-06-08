
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
function (op::Operator{:∫vᵢσdΩ_mohr_coulomb})(ap::T,k::AbstractMatrix{Float64},fint::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    𝜙 = op.𝜙
    tol = op.tol 
    Cᵢᵢᵢᵢ = λ + 2*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    C₁₁₁₁ = Cᵢᵢᵢᵢ
    C₂₂₂₂ = Cᵢᵢᵢᵢ
    C₃₃₃₃ = Cᵢᵢᵢᵢ
    C₁₁₂₂ = Cᵢᵢⱼⱼ
    C₁₁₃₃ = Cᵢᵢⱼⱼ
    C₂₂₃₃ = Cᵢᵢⱼⱼ
    C₁₂₁₂ = Cᵢⱼᵢⱼ
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₃₃ = ξ.σ₃₃
        σ₁₂ = ξ.σ₁₂
        sₙ = ξ.sₙ
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
        σ₁₁ᵗʳ =σ₁₁+C₁₁₁₁*Δε₁₁+C₁₂₁₂*Δε₂₂
        σ₂₂ᵗʳ =σ₂₂+C₂₂₂₂*Δε₂₂+C₁₂₁₂*Δε₁₁
        σ₃₃ᵗʳ =0
        σ₁₂ᵗʳ =σ₁₂+2.0*C₁₂₁₂*Δε₁₂ 
        σ₁,σ₂,σ₃,n₁,n₂,n₃ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₃₃ᵗʳ,σ₁₂ᵗʳ,0.0,0.0)
        fᵗʳ = σ₁-σ₃ + (σ₁+σ₃) * sin(𝜙) - 2*c*cos(𝜙)
        if fᵗʳ > tol
            Δγ = fᵗʳ/(4.0*sin(𝜙)^2*λ+4.0*(1.0+sin(𝜙)^2)*μ)
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*(C₁₁₁₁*((sin𝜙+1)*n₁⁽¹⁾*n₁⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₁⁽³⁾)+C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*(C₂₂₂₂*((sin𝜙+1)*n₂⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₂⁽³⁾*n₂⁽³⁾)+C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*2.0*(C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*(C₁₁₁₁*((sin𝜙+1)*n₁⁽¹⁾*n₁⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₁⁽³⁾)+C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*(C₂₂₂₂*((sin𝜙+1)*n₂⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₂⁽³⁾*n₂⁽³⁾)+C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*2.0*(C₁₂₁₂*((sin𝜙+1)*n₁⁽¹⁾*n₂⁽¹⁾+(sin𝜙-1)n₁⁽³⁾*n₂⁽³⁾))
            sin²𝜙 = sin(𝜙)^2
            sin²𝜙psin𝜙 = sin(𝜙)^2+sin(𝜙) 
            sin²𝜙dsin𝜙=sin(𝜙)^2-sin(𝜙) 
            Cᵗ₁₁₁₁=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*n₁⁽¹⁾^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*n₁⁽³⁾^2-8.0cos(𝜙)^2*μ^2*n₁⁽¹⁾^2*n₁⁽³⁾^2+4.0*(sin(𝜙)+1)^2*μ^2*n₁⁽¹⁾^4+4.0*(sin(𝜙)-1)^2*μ^2*n₁⁽³⁾^4)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗ₂₂₂₂=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*n₂⁽¹⁾^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*n₂⁽³⁾^2-8.0cos(𝜙)^2*μ^2*n₂⁽¹⁾^2*n₂⁽³⁾^2+4.0*(sin(𝜙)+1)^2*μ^2*n₂⁽¹⁾^4+4.0*(sin(𝜙)-1)^2*μ^2*n₂⁽³⁾^4)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗ₃₃₃₃=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*n₃⁽¹⁾^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*n₃⁽³⁾^2-8.0cos(𝜙)^2*μ^2*n₃⁽¹⁾^2*n₃⁽³⁾^2+4.0*(sin(𝜙)+1)^2*μ^2*n₃⁽¹⁾^4+4.0*(sin(𝜙)-1)^2*μ^2*n₃⁽³⁾^4)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗ₁₁₂₂=Cᵢᵢⱼⱼ-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*(n₁⁽¹⁾^2+n₂⁽¹⁾^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*(n₁⁽³⁾^2+n₂⁽³⁾^2)+4.0*(sin(𝜙)+1)^2*μ^2*n₁⁽¹⁾^2*n₂⁽¹⁾^2-4.0cos(𝜙)^2*μ^2*(n₁⁽³⁾^2*n₂⁽¹⁾^2+n₁⁽¹⁾^2*n₂⁽³⁾^2)+4.0*(sin(𝜙)-1)^2*μ^2*n₁⁽³⁾^2*n₂⁽³⁾^2)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗ₁₁₃₃=Cᵢᵢⱼⱼ-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*(n₁⁽¹⁾^2+n₃⁽¹⁾^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*(n₁⁽³⁾^2+n₃⁽³⁾^2)+4.0*(sin(𝜙)+1)^2*μ^2*n₁⁽¹⁾^2*n₃⁽¹⁾^2-4.0cos(𝜙)^2*μ^2*(n₁⁽³⁾^2*n₃⁽¹⁾^2+n₁⁽¹⁾^2*n₃⁽³⁾^2)+4.0*(sin(𝜙)-1)^2*μ^2*n₁⁽³⁾^2*n₃⁽³⁾^2)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗ₂₂₃₃=Cᵢᵢⱼⱼ-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*(n₂⁽¹⁾^2+n₂⁽¹⁾^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂⁽³⁾^2+n₂⁽³⁾^2)+4.0*(sin(𝜙)+1)^2*μ^2*n₂⁽¹⁾^2*n₂⁽¹⁾^2-4.0cos(𝜙)^2*μ^2*(n₂⁽³⁾^2*n₂⁽¹⁾^2+n₂⁽¹⁾^2*n₂⁽³⁾^2)+4.0*(sin(𝜙)-1)^2*μ^2*n₂⁽³⁾^2*n₂⁽³⁾^2)/4(sin(𝜙)^2*λ+4(1.0+sin(𝜙)^2)*μ)
            Cᵗᵢⱼᵢⱼ = Cᵢⱼᵢⱼ
        else
            ξ.σₙ = σᵗʳ
            Cᵗᵢᵢᵢᵢ = Cᵢᵢᵢᵢ
            Cᵗᵢᵢⱼⱼ = Cᵢᵢⱼⱼ
            Cᵗᵢⱼᵢⱼ = Cᵢⱼᵢⱼ           

        end
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J-1] += (Cᵗ₂₂₂₂*B₁[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J-1] += (Cᵗ₃₃₃₃*B₁[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤 
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₃₃*B₁[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I-1,2*J]   += (Cᵗ₂₂₃₃*B₁[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₃₃*B₂[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₂₂₃₃*B₂[i]*B₁[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₁₁₁₁*B₂[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₃₃₃₃*B₂[i]*B₂[j] + Cᵗᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
            f[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            f[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
        end
    end
end    