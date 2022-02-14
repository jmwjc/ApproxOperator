
function checkIC(a::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},b::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = a.𝓒[1].x;y₁ = a.𝓒[1].y
    x₂ = a.𝓒[2].x;y₂ = a.𝓒[2].y
    x₃ = a.𝓒[3].x;y₃ = a.𝓒[3].y
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    nᵈ = length(a.𝓒)
    # nᵖ = get𝑛𝒒(a)
    nᵖ = get𝑛𝒑(a)
    f₁ = zeros(nᵈ,nᵖ)
    f₂ = zeros(nᵈ,nᵖ)
    for ξ in a.𝓖
        N,B₁,B₂ = get∇𝝭(a,ξ)
        𝒒 = get𝒒(a,ξ)
        p = get𝒑(a,ξ)
        𝑤 = get𝑤(a,ξ)
        for i in 1:nᵈ
            for j in 1:nᵖ
                # f₁[i,j] += B₁[i]*𝒒[j]*𝑤
                # f₂[i,j] += B₂[i]*𝒒[j]*𝑤
                f₁[i,j] += B₁[i]*p[j]*𝑤
                f₂[i,j] += B₂[i]*p[j]*𝑤
            end
        end
    end
    for ξ in b.𝓖
        N = get𝝭(b,ξ)
        𝒒,∂𝒒∂ξ,∂𝒒∂η = get∇𝒒(b,ξ)
        𝒙 = get𝒙(b,ξ)
        p,∂p∂x,∂p∂y = get∇𝒑(b,𝒙)
        wᵇ = ξ.wᵇ
        w = ξ.w
        𝑤 = get𝑤(b,ξ)
        nᵇ₁ = 0.0;nᵇ₂ = 0.0
        nᵇ₁ += ξ.ξ == 0.0 ? n₁₁ : 0.0
        nᵇ₁ += ξ.η == 0.0 ? n₂₁ : 0.0
        nᵇ₁ += ξ.ξ+ξ.η ≈ 1.0 ? n₃₁ : 0.0
        nᵇ₂ += ξ.ξ == 0.0 ? n₁₂ : 0.0
        nᵇ₂ += ξ.η == 0.0 ? n₂₂ : 0.0
        nᵇ₂ += ξ.ξ+ξ.η ≈ 1.0 ? n₃₂ : 0.0
        for j in 1:nᵖ
            # W₁ = 𝒒[j]*nᵇ₁*wᵇ + (∂𝒒∂ξ[j]*n₁₁+∂𝒒∂η[j]*n₂₁)*w
            # W₂ = 𝒒[j]*nᵇ₂*wᵇ + (∂𝒒∂ξ[j]*n₁₂+∂𝒒∂η[j]*n₂₂)*w
            W₁ = p[j]*nᵇ₁*wᵇ + ∂p∂x[j]*𝑤
            W₂ = p[j]*nᵇ₂*wᵇ + ∂p∂y[j]*𝑤
            for i in 1:nᵈ
                f₁[i,j] -= N[i]*W₁
                f₂[i,j] -= N[i]*W₂
            end
        end
    end
    return f₁,f₂
end

function checkCC(a::ReproducingKernel)
    nᵖ = get𝑛𝒑(a)
    nⁱ = length(a.𝓖)
    f = zeros(nⁱ,nᵖ)
    for i in 1:length(a.𝓖)
        ξ = a.𝓖[i]
        𝒙 = get𝒙(a,ξ)
        p = get𝒑(a,𝒙)
        N,B₁,B₂ = get∇𝝭(a,ξ)
        for j in 1:nᵖ
            for k in 1:length(a.𝓒)
                𝒙̄ = (a.𝓒[k].x,a.𝓒[k].y,a.𝓒[k].z)
                p̄ = get𝒑(a,𝒙̄)
                f[i,j] += N[k]*p̄[j]
            end
            f[i,j] -= p[j]
        end
    end
    return f
end
