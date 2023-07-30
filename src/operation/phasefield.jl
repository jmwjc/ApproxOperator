
function (op::Operator{:∫v²uₓuₓdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*B[i]*B[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx_hard_device})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
        end
        ℋₜ = max(ℋ,(ε-1)^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*kc/2/l*𝑤
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        v = 0.0
        ∂v∂x = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
            ∂v∂x += B[i]*xᵢ.v
        end
        ℋₜ = max(ℋ,(v^2+η)*ε^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*(kc/2/l - η*ℋₜ)*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_1D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,(v^2+η)*ε^2)
    end
end