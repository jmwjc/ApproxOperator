
function (op::Operator{:L₂})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        ūᵢ = ξ.u
        uᵢ = 0
        for (i,xᵢ) in enumerate(ap.𝓒)
            uᵢ += N[i]*xᵢ.d
        end
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū²  += ūᵢ^2*𝑤
    end
    return Δu², ū²
end

function (op::Operator{:L₂})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:H₁})(ap::T) where T<:AbstractElement
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂ūᵢ∂z = ξ.∂u∂z
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂uᵢ∂z = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂uᵢ∂z += B₃[i]*xᵢ.d
        end
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2 + (∂uᵢ∂z - ∂ūᵢ∂z)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2 + ∂ūᵢ∂z^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇u², ∇ū², Δu², ū²
end

function (op::Operator{:H₁})(aps::Vector{T}) where T<:AbstractElement
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:Hₑ_PlaneStress})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂
        σ̄₂₂ = Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₁₁
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_PlaneStress})(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:Hₑ_Incompressible})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵈ = E/(1+ν)
    for ξ in ap.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = 0.5*(∂ū₁∂y + ∂ū₂∂x)
        σ̄₁₁ = Cᵈ*(2*ε̄₁₁ -   ε̄₂₂)/3
        σ̄₂₂ = Cᵈ*(- ε̄₁₁ + 2*ε̄₂₂)/3
        σ̄₁₂ = Cᵈ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += 0.5*(B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂)
        end
        σ₁₁ = Cᵈ*(2*ε₁₁ -   ε₂₂)/3
        σ₂₂ = Cᵈ*(- ε₁₁ + 2*ε₂₂)/3
        σ₁₂ = Cᵈ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_Incompressible})(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function set∇𝑢!(ap::T) where T<:AbstractElement
    for ξ in ap.𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x_]
        B₂ = ξ[:∂𝝭∂y_]
        B₃ = ξ[:∂𝝭∂z_]
        𝒙 = (ξ.x,ξ.y,ξ.z)
        u = 0.
        ∂u∂x = 0.
        ∂u∂y = 0.
        ∂u∂z = 0.
        for (i,xᵢ) in enumerate(ap.𝓒)
            u += N[i]*x.d
            ∂u∂x += B₁[i]*x.d
            ∂u∂y += B₂[i]*x.d
            ∂u∂z += B₃[i]*x.d
        end
        ξ.x = 𝒙[1]
        ξ.y = 𝒙[2]
        ξ.z = 𝒙[3]
        ξ.u = u
        ξ.∂u∂x = ∂u∂x
        ξ.∂u∂y = ∂u∂y
        ξ.∂u∂z = ∂u∂z
    end
end

function (op::Operator{:H₃})(aps::Vector{T}) where T<:AbstractElement
    H₃Norm_Δu²= 0
    H₃Norm_ū² = 0
    H₂Norm_Δu²= 0
    H₂Norm_ū² = 0
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇³u², ∇³ū²,Δ∇²u², ∇²ū²,Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₃Norm_Δu² += Δu² + Δ∇u² + Δ∇²u² + Δ∇³u²
        H₃Norm_ū²  += ū² + ∇ū² + ∇²ū² + ∇³ū²
        H₂Norm_Δu² += Δu² + Δ∇u² + Δ∇²u²
        H₂Norm_ū²  += ū² + ∇ū² + ∇²ū²
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₃Norm_Δu²/H₃Norm_ū²)^0.5, (H₂Norm_Δu²/H₂Norm_ū²)^0.5, (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function (op::Operator{:H₃})(ap::T) where T<:AbstractElement
    Δ∇³u²= 0
    ∇³ū² = 0
    Δ∇²u²= 0
    ∇²ū² = 0
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B₁₁₁ = ξ[:∂³𝝭∂x³]
        B₁₁₂ = ξ[:∂³𝝭∂x²∂y]
        B₁₂₂ = ξ[:∂³𝝭∂x∂y²]
        B₂₂₂ = ξ[:∂³𝝭∂y³]
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂²ūᵢ∂x² = ξ.∂²u∂x²
        ∂²ūᵢ∂x∂y = ξ.∂²u∂x∂y
        ∂²ūᵢ∂y² = ξ.∂²u∂y²
        ∂³ūᵢ∂x³ = ξ.∂³u∂x³
        ∂³ūᵢ∂x²∂y = ξ.∂³u∂x²∂y
        ∂³ūᵢ∂x∂y² = ξ.∂³u∂x∂y²
        ∂³ūᵢ∂y³ = ξ.∂³u∂y³
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂²uᵢ∂x² = 0.
        ∂²uᵢ∂x∂y = 0.
        ∂²uᵢ∂y² = 0.
        ∂³uᵢ∂x³ = 0.
        ∂³uᵢ∂x²∂y = 0.
        ∂³uᵢ∂x∂y² = 0.
        ∂³uᵢ∂y³ = 0.
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂²uᵢ∂x² += B₁₁[i]*xᵢ.d
            ∂²uᵢ∂x∂y += B₁₂[i]*xᵢ.d
            ∂²uᵢ∂y² += B₂₂[i]*xᵢ.d
            ∂³uᵢ∂x³ += B₁₁₁[i]*xᵢ.d
            ∂³uᵢ∂x²∂y += B₁₁₂[i]*xᵢ.d
            ∂³uᵢ∂x∂y² += B₁₂₂[i]*xᵢ.d
            ∂³uᵢ∂y³ += B₂₂₂[i]*xᵢ.d
        end
        Δ∇³u² += ((∂³uᵢ∂x³ - ∂³ūᵢ∂x³)^2 + (∂³uᵢ∂x²∂y - ∂³ūᵢ∂x²∂y)^2 + (∂³uᵢ∂x∂y² - ∂³ūᵢ∂x∂y²)^2 + (∂³uᵢ∂y³ - ∂³ūᵢ∂y³)^2)*𝑤
        ∇³ū² += (∂³ūᵢ∂x³^2 + ∂³ūᵢ∂x²∂y^2  + ∂³ūᵢ∂x∂y²^2+ ∂³ūᵢ∂y³^2)*𝑤
        Δ∇²u² += ((∂²uᵢ∂x² - ∂²ūᵢ∂x²)^2 + (∂²uᵢ∂x∂y - ∂²ūᵢ∂x∂y)^2 + (∂²uᵢ∂y² - ∂²ūᵢ∂y²)^2)*𝑤
        ∇²ū² += (∂²ūᵢ∂x²^2 + ∂²ūᵢ∂x∂y^2 + ∂²ūᵢ∂y²^2)*𝑤
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇³u², ∇³ū², Δ∇²u², ∇²ū², Δ∇u², ∇ū², Δu², ū²
end