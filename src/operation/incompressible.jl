
function (op::Operator{:∫∫qpdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    K = E/3/(1-2*ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] -= N[i]*N[j]/K*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫qpdxdy})(aₚ::T,aₛ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ₚ = aₚ.𝓒
    𝓒ₛ = aₛ.𝓒
    𝓖ₚ = aₚ.𝓖
    𝓖ₛ = aₛ.𝓖
    E = op.E
    ν = op.ν
    K = E/3/(1-2*ν)
    for (ξₚ,ξₛ) in zip(𝓖ₚ,𝓖ₛ)
        𝑤 = ξₚ.𝑤
        Nₚ = ξₚ[:𝝭]
        Nₛ = ξₛ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ₚ)
                J = xⱼ.𝐼
                k[I,J] -= N[i]*N[j]/K*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫p∇vdxdy})(aᵤ::T,aₚ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] += N[i]*B₁[j]*𝑤
                k[I,2*J]   += N[i]*B₂[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫pnᵢgᵢds})(aᵤ::T,aₚ::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        Nₚ = ξₚ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[I,2*J-1] -= Nₚ[i]*Nᵤ[j]*(n₁*n₁₁+n₂*n₁₂)*𝑤
                k[I,2*J]   -= Nₚ[i]*Nᵤ[j]*(n₁*n₁₂+n₂*n₂₂)*𝑤
            end
            f[I] -= Nₚ[i]*(n₁*n₁₁*g₁+n₁*n₁₂*g₂+n₂*n₁₂*g₁+n₂*n₂₂*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫∫δsᵢⱼsᵢⱼdxdy})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    G = E/(1+ν)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[4*I-3,4*J-3] -= 4/9/G*N[i]*N[j]*𝑤
                k[4*I-3,4*J-2] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-3,4*J-1] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-2,4*J-3] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-2,4*J-2] -= 2/3/G*N[i]*N[j]*𝑤
                k[4*I-2,4*J-1] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-1,4*J-3] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-1,4*J-2] += 1/3/G*N[i]*N[j]*𝑤
                k[4*I-1,4*J-1] -= 2/3/G*N[i]*N[j]*𝑤
                k[4*I,4*J]     -=   2/G*N[i]*N[j]*𝑤

                # k[3*I-2,3*J-2] -= N[i]*N[j]/G*𝑤
                # k[3*I-1,3*J-1] -= N[i]*N[j]/G*𝑤
                # k[3*I,3*J]     -= 2*N[i]*N[j]/G*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫sᵢⱼεᵢⱼdxdy})(aᵤ::T,aₛ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₛ = aₛ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        N = ξₛ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] += 2/3*N[i]*B₁[j]*𝑤
                k[4*I-3,2*J]   -= 1/3*N[i]*B₂[j]*𝑤
                k[4*I-2,2*J-1] -= 1/3*N[i]*B₁[j]*𝑤
                k[4*I-2,2*J]   += 2/3*N[i]*B₂[j]*𝑤
                k[4*I-1,2*J-1] -= 1/3*N[i]*B₁[j]*𝑤
                k[4*I-1,2*J]   -= 1/3*N[i]*B₂[j]*𝑤
                k[4*I,2*J-1]   +=     N[i]*B₂[j]*𝑤
                k[4*I,2*J]     +=     N[i]*B₁[j]*𝑤

                # k[3*I-2,2*J-1] += 2/3*N[i]*B₁[j]*𝑤
                # k[3*I-2,2*J]   -= 1/3*N[i]*B₂[j]*𝑤
                # k[3*I-1,2*J-1] -= 1/3*N[i]*B₁[j]*𝑤
                # k[3*I-1,2*J]   += 2/3*N[i]*B₂[j]*𝑤
                # k[3*I,2*J-1]   += N[i]*B₂[j]*𝑤
                # k[3*I,2*J]     += N[i]*B₁[j]*𝑤

                # k[3*I-2,2*J-1] += N[i]*B₁[j]*𝑤
                # k[3*I-1,2*J]   += N[i]*B₂[j]*𝑤
                # k[3*I,2*J-1]   += N[i]*B₂[j]*𝑤
                # k[3*I,2*J]     += N[i]*B₁[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫sᵢⱼnⱼgᵢds})(aᵤ::T,aₛ::S;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₛ = aₛ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₛ = aₛ.𝓖
    for (ξᵤ,ξₛ) in zip(𝓖ᵤ,𝓖ₛ)
        Nₛ = ξₛ[:𝝭]
        Nᵤ = ξᵤ[:𝝭]
        n₁ = ξᵤ.n₁
        n₂ = ξᵤ.n₂
        g₁ = ξᵤ.g₁
        g₂ = ξᵤ.g₂
        n₁₁ = ξᵤ.n₁₁
        n₁₂ = ξᵤ.n₁₂
        n₂₂ = ξᵤ.n₂₂
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₛ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[4*I-3,2*J-1] -= Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₁ - 1/3*n₂*n₁₂)*𝑤
                k[4*I-3,2*J]   -= Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₂ - 1/3*n₂*n₂₂)*𝑤
                k[4*I-2,2*J-1] -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₁ + 2/3*n₂*n₁₂)*𝑤
                k[4*I-2,2*J]   -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₂ + 2/3*n₂*n₂₂)*𝑤
                k[4*I-1,2*J-1] -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₁ - 1/3*n₂*n₁₂)*𝑤
                k[4*I-1,2*J]   -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₂ - 1/3*n₂*n₂₂)*𝑤
                k[4*I,2*J-1]   -= Nₛ[i]*Nᵤ[j]*(n₁*n₁₂ + n₂*n₁₁)*𝑤
                k[4*I,2*J]     -= Nₛ[i]*Nᵤ[j]*(n₁*n₂₂ + n₂*n₁₂)*𝑤

                # k[3*I-2,2*J-1] -= Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₁ - 1/3*n₂*n₁₂)*𝑤
                # k[3*I-2,2*J]   -= Nₛ[i]*Nᵤ[j]*( 2/3*n₁*n₁₂ - 1/3*n₂*n₂₂)*𝑤
                # k[3*I-1,2*J-1] -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₁ + 2/3*n₂*n₁₂)*𝑤
                # k[3*I-1,2*J]   -= Nₛ[i]*Nᵤ[j]*(-1/3*n₁*n₁₂ + 2/3*n₂*n₂₂)*𝑤
                # k[3*I,2*J-1]   -= Nₛ[i]*Nᵤ[j]*(n₁*n₁₂+n₂*n₁₁)*𝑤
                # k[3*I,2*J]     -= Nₛ[i]*Nᵤ[j]*(n₁*n₂₂+n₂*n₁₂)*𝑤
            end
            f[4*I-3] -= Nₛ[i]*(( 2/3*n₁*n₁₁-1/3*n₂*n₁₂)*g₁+( 2/3*n₁*n₁₂-1/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I-2] -= Nₛ[i]*((-1/3*n₁*n₁₁+2/3*n₂*n₁₂)*g₁+(-1/3*n₁*n₁₂+2/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I-1] -= Nₛ[i]*((-1/3*n₁*n₁₁-1/3*n₂*n₁₂)*g₁+(-1/3*n₁*n₁₂-1/3*n₂*n₂₂)*g₂)*𝑤
            f[4*I]   -= Nₛ[i]*((n₁*n₁₂+n₂*n₁₁)*g₁+(n₁*n₂₂+n₂*n₁₂)*g₂)*𝑤

            # f[3*I-2] -= Nₛ[i]*((2/3*n₁*n₁₁-1/3*n₂*n₁₂)*g₁+(2/3*n₁*n₁₂-1/3*n₂*n₂₂)*g₂)*𝑤
            # f[3*I-1] -= Nₛ[i]*((-1/3*n₁*n₁₁+2/3*n₂*n₁₂)*g₁+(-1/3*n₁*n₁₂+2/3*n₂*n₂₂)*g₂)*𝑤
            # f[3*I]   -= Nₛ[i]*((n₁*n₁₂+n₂*n₁₁)*g₁+(n₁*n₂₂+n₂*n₁₂)*g₂)*𝑤
        end
    end
end


