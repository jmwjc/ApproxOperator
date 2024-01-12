

function (op::Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    h = op.h
    Dᵐ = E*h/(1-ν^2)
    Dᵇ = E*h^3/12/(1-ν^2)
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        a¹¹ = ξ.a¹¹
        a¹² = ξ.a¹²
        a²² = ξ.a²²
        a₁₁ = ξ.a₁₁
        a₁₂ = ξ.a₁₂
        a₁₃ = ξ.a₁₃
        a₂₁ = ξ.a₂₁
        a₂₂ = ξ.a₂₂
        a₂₃ = ξ.a₂₃
        a₃₁ = ξ.a₃₁
        a₃₂ = ξ.a₃₂
        a₃₃ = ξ.a₃₃
        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₂₂ = ξ.Γ²₂₂
        Γ²₁₂ = ξ.Γ²₁₂
        𝐃 = @SArray [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
                     ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
                                       a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ᵐᵢ = @SArray [          a₁₁*B₁[i]           a₁₂*B₁[i]           a₁₃*B₁[i];
                                     a₂₁*B₂[i]           a₂₂*B₂[i]           a₂₃*B₂[i];
                           a₁₁*B₂[i]+a₂₁*B₁[i] a₁₂*B₂[i]+a₂₂*B₁[i] a₁₃*B₂[i]+a₂₃*B₁[i]]
            Bᵇ₁₁ = Γ¹₁₁*B₁[i]+Γ²₁₁*B₂[i]-B₁₁[i]
            Bᵇ₂₂ = Γ¹₂₂*B₁[i]+Γ²₂₂*B₂[i]-B₂₂[i]
            Bᵇ₁₂ = Γ¹₁₂*B₁[i]+Γ²₁₂*B₂[i]-B₁₂[i]
            𝐁ᵇᵢ = @SArray [  Bᵇ₁₁*a₃₁   Bᵇ₁₁*a₃₂   Bᵇ₁₁*a₃₃;
                             Bᵇ₂₂*a₃₁   Bᵇ₂₂*a₃₂   Bᵇ₂₂*a₃₃;
                           2*Bᵇ₁₂*a₃₁ 2*Bᵇ₁₂*a₃₂ 2*Bᵇ₁₂*a₃₃]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ᵐⱼ = @SArray [          a₁₁*B₁[j]           a₁₂*B₁[j]           a₁₃*B₁[j];
                                         a₂₁*B₂[j]           a₂₂*B₂[j]           a₂₃*B₂[j];
                               a₁₁*B₂[j]+a₂₁*B₁[j] a₁₂*B₂[j]+a₂₂*B₁[j] a₁₃*B₂[j]+a₂₃*B₁[j]]
                Bᵇ₁₁ = Γ¹₁₁*B₁[j]+Γ²₁₁*B₂[j]-B₁₁[j]
                Bᵇ₂₂ = Γ¹₂₂*B₁[j]+Γ²₂₂*B₂[j]-B₂₂[j]
                Bᵇ₁₂ = Γ¹₁₂*B₁[j]+Γ²₁₂*B₂[j]-B₁₂[j]
                𝐁ᵇⱼ = @SArray [  Bᵇ₁₁*a₃₁   Bᵇ₁₁*a₃₂   Bᵇ₁₁*a₃₃;
                                 Bᵇ₂₂*a₃₁   Bᵇ₂₂*a₃₂   Bᵇ₂₂*a₃₃;
                               2*Bᵇ₁₂*a₃₁ 2*Bᵇ₁₂*a₃₂ 2*Bᵇ₁₂*a₃₃]
                kᵐ = Dᵐ*𝐁ᵐᵢ'*𝐃*𝐁ᵐⱼ
                kᵇ = Dᵇ*𝐁ᵇᵢ'*𝐃*𝐁ᵇⱼ
                k[3*I-2,3*J-2] += (kᵐ[1,1] + kᵇ[1,1])*𝑤
                k[3*I-2,3*J-1] += (kᵐ[1,2] + kᵇ[1,2])*𝑤
                k[3*I-2,3*J]   += (kᵐ[1,3] + kᵇ[1,3])*𝑤
                k[3*I-1,3*J-2] += (kᵐ[2,1] + kᵇ[2,1])*𝑤
                k[3*I-1,3*J-1] += (kᵐ[2,2] + kᵇ[2,2])*𝑤
                k[3*I-1,3*J]   += (kᵐ[2,3] + kᵇ[2,3])*𝑤
                k[3*I,3*J-2]   += (kᵐ[3,1] + kᵇ[3,1])*𝑤
                k[3*I,3*J-1]   += (kᵐ[3,2] + kᵇ[3,2])*𝑤
                k[3*I,3*J]     += (kᵐ[3,3] + kᵇ[3,3])*𝑤
            end
        end
    end
end

function (op::Operator{:∫δθθdΓ})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n¹ = ξ.n₁
        n² = ξ.n₂
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            𝐁ₙᵢ = @SArray [(B₁[i]*n¹ + B₂[i]*n²)*a₃₁;(B₁[i]*n¹ + B₂[i]*n²)*a₃₂;(B₁[i]*n¹ + B₂[i]*n²)*a₃₃]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                𝐁ₙⱼ = @SArray [(B₁[j]*n¹ + B₂[j]*n²)*a₃₁;(B₁[j]*n¹ + B₂[j]*n²)*a₃₂;(B₁[j]*n¹ + B₂[j]*n²)*a₃₃]
                k[3*I-2,3*J-2] += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[1]*𝑤
                k[3*I-2,3*J-1] += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[2]*𝑤
                k[3*I-2,3*J]   += α*𝐁ₙᵢ[1]*𝐁ₙⱼ[3]*𝑤
                k[3*I-1,3*J-2] += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[1]*𝑤
                k[3*I-1,3*J-1] += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[2]*𝑤
                k[3*I-1,3*J]   += α*𝐁ₙᵢ[2]*𝐁ₙⱼ[3]*𝑤
                k[3*I,3*J-2]   += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[1]*𝑤
                k[3*I,3*J-1]   += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[2]*𝑤
                k[3*I,3*J]     += α*𝐁ₙᵢ[3]*𝐁ₙⱼ[3]*𝑤
            end
            f[3*I-2] += α*𝐁ₙᵢ[1]*θ*𝑤
            f[3*I-1] += α*𝐁ₙᵢ[2]*θ*𝑤
            f[3*I]   += α*𝐁ₙᵢ[3]*θ*𝑤
        end
    end
end

function (op::Operator{:ScordelisLoRoof_𝐴})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒
    ξ = ap.𝓖[1]
    w = 0.0
    N = ξ[:𝝭]
    for (i,xᵢ) in enumerate(𝓒)
        w += N[i]*xᵢ.d₃
    end
    return w
end
