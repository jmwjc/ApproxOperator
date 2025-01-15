module Hamilton
    
using ..ApproxOperator: AbstractElement

function ∫∫q̇mṗqkpdxdt(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        m = ξ.m
        kᶜ = ξ.k
        B = ξ[:∂𝝭∂x]
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (B[i]*m*B[j] - N[i]*kᶜ*N[j])*𝑤
            end
        end
    end
end

function ∫∫∇q∇pdxdt(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        EA = ξ.EA
        Bₓ = ξ[:∂𝝭∂x]
        Bₜ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (-Bₜ[i]*ρA*Bₜ[j] + Bₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

function ∫∫∇q∇pdxdt(a₁::T,a₂::S,k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒₁ = a₁.𝓒; 𝓖₁ = a₁.𝓖
    𝓒₂ = a₂.𝓒; 𝓖₂ = a₂.𝓖
    for (ξ₁,ξ₂) in (𝓖₁,𝓖₂)
        B̄ₓ = ξ₁[:∂𝝭∂x]
        B̄ₜ = ξ₁[:∂𝝭∂y]
        ρA = ξ₂.ρA
        EA = ξ₂.EA
        Bₓ = ξ₂[:∂𝝭∂x]
        Bₜ = ξ₂[:∂𝝭∂y]
        𝑤 = ξ₂.𝑤
        for (i,xᵢ) in enumerate(𝓒₁)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒₂)
                J = xⱼ.𝐼
                k[I,J] += (-B̄ₜ[i]*ρA*Bₜ[j] + B̄ₓ[i]*EA*Bₓ[j])*𝑤
            end
        end
    end
end

function stabilization_bar_LSG(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ρA = ξ.ρA
        EA = ξ.EA
        Β = ξ.Β
        Bₓₓ = ξ[:∂²𝝭∂x²]
        Bₜₜ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += Β*(ρA*Bₜₜ[i] - EA*Bₓₓ[i])*(ρA*Bₜₜ[j] - EA*Bₓₓ[j])*𝑤
            end
        end
    end
end

end