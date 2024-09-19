function (op::Operator{:∫∫pᵢD⁻¹pⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += N[i]*N[j]*𝑤
                k[2*I-1,2*J]   += 0
                k[2*I,2*J-1]   += 0
                k[2*I,2*J]     += N[i]*N[j]*𝑤
            end 
        end
    end
end
function (op::Operator{:∫∫∇TᵢD∇Tⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    t = op.t
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (D*t*B₁[i]*B₁[j]+D*t*B₂[i]*B₂[j])*𝑤 
            end 
        end
    end
end
function (op::Operator{:∫∫pᵢ∇uⱼdxdy})(aₚ::T,aᵤ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ₚ)
                J = xⱼ.𝐼
                k[I,2*J-1] += -N[j]*B₁[i]*𝑤
                k[I,2*J]   += -N[j]*B₂[i]*𝑤
            end
        end
    end
end
function (op::Operator{:∫∫∇𝒑udxdy})(aₚ::T,aᵤ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ  = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξᵤ[:𝝭]
        B₁ = ξₚ[:∂𝝭∂x]
        B₂ = ξₚ[:∂𝝭∂y]
        𝑤 = ξₚ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,J] += B₁[i]*N[j]*𝑤
                k[2*I,J]   += B₂[i]*N[j]*𝑤
            end
        end
    end
end
function (op::Operator{:∫pᵢnᵢuds})(aₚ::T,aᵤ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒 
    𝓒ₚ = aₚ.𝓒  
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
  
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        Nᵤ = ξᵤ[:𝝭]
        Nₚ = ξₚ[:𝝭]
        𝑤 = ξᵤ.𝑤
        n₁ = ξₚ.n₁
        n₂ = ξₚ.n₂
        g = ξᵤ.g
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,J] -= Nₚ[i]*Nᵤ[j]*n₁*𝑤
                k[2*I,J]   -= Nₚ[i]*Nᵤ[j]*n₂*𝑤
            end
        end
    end
end

function (op::Operator{:∫pᵢnᵢgⱼds})(aₚ::T,aᵤ::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒 
    𝓒ₚ = aₚ.𝓒  
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
  
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        Nᵤ = ξᵤ[:𝝭]
        Nₚ = ξₚ[:𝝭]
        𝑤 = ξᵤ.𝑤
        n₁ = ξₚ.n₁
        n₂ = ξₚ.n₂
        g = ξᵤ.g
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,J] += Nₚ[i]*Nᵤ[j]*n₁*𝑤
                k[2*I,J]   += Nₚ[i]*Nᵤ[j]*n₂*𝑤
            end
            f[2*I-1] += Nₚ[i]*n₁*g*𝑤
            f[2*I] += Nₚ[i]*n₂*g*𝑤
        end
    end
end

function (op::Operator{:∫∫Tᵢsᵢdxdy})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    t = op.t
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        s = ξ.s
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += t*N[i]*s*𝑤
        end
    end
end
function (op::Operator{:∫Tᵢhᵢds})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        h = ξ.h
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*h*𝑤
        end
    end
end

function (op::Operator{:∫Tᵢgᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    t = op.t
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*t*N[i]*N[j]*𝑤
            end
            f[I] += α*t*N[i]*g*𝑤
        end
    end
end