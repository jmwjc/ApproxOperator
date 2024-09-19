
function (op::Operator{:∫ESdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖 
    Eᵉ = op.E
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        F = 1.0
        for (i,xᵢ) in enumerate(𝓒)
            F += B[i]*xᵢ.d
        end
        E = 0.5*(F^2-1.0)
        S = Eᵉ*E
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*(S+F^2*Eᵉ)*B[j]*𝑤
            end
            f[I] += B[i]*Eᵉ*F*E*𝑤
        end
    end
end
function (op::Operator{:∫∫pJdxdy})(aᵤ::T,aₚ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        𝑤 = ξᵤ.𝑤
        for (i,xᵢ) in  enumerate(𝓒ᵤ)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        # J = F₁₁*F₂₂-F₁₂*F₂₁
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # k[I,2*J-1] += J*N[i]*B₁[j]*𝑤
                # k[I,2*J]   += J*N[i]*B₂[j]*𝑤
                k[I,2*J-1] += N[i]*B₁[j]*𝑤
                k[I,2*J]   += N[i]*B₂[j]*𝑤
            end
        end
    end
end
function (op::Operator{:∫∫pCdxy_incompressible})(aᵤ::T,aₚ::S;f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
 
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        N = ξₚ[:𝝭]
        𝑤 = ξᵤ.𝑤
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        J = F₁₁*F₂₂-F₁₂*F₂₁
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼
            f[I]-=(N[i]*(J-1.0))*𝑤
        end
    end
end

function (op::Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        Cᵢⱼᵢⱼ = E/(1+ν)/2
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        
        J = F₁₁*F₂₂-F₁₂*F₂₁
        F₁₁ = J^(-1/3)*F₁₁
        F₂₂ = J^(-1/3)*F₂₂
        F₁₂ = J^(-1/3)*F₁₂
        F₂₁ = J^(-1/3)*F₂₁

        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +  (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ
                               +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤
                              
                k[2*I-1,2*J]   += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ)*𝑤
                               
                k[2*I,2*J-1]   += ((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ)*𝑤

                k[2*I,2*J]     += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ 
                               +   (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ
                               +    B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂)*𝑤

            end
        end
    end
end
function (op::Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff_mix})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        # Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        # Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        # Cᵢⱼᵢⱼ = E/(1+ν)/2

        λ = E*ν/(1+ν)/(1-2*ν)
        μ = E/(1+ν)/2
        Cᵢᵢᵢᵢ = 2*μ
        Cᵢⱼᵢⱼ = μ
        Cᵢᵢⱼⱼ = 0
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        
        J = F₁₁*F₂₂-F₁₂*F₂₁
        # F₁₁ = J^(-1/3)*F₁₁
        # F₂₂ = J^(-1/3)*F₂₂
        # F₁₂ = J^(-1/3)*F₁₂
        # F₂₁ = J^(-1/3)*F₂₁


        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁
        S₂₂ = Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        I₁ = F₁₁*F₁₁+F₂₂*F₂₂
        I₂ = I₁^2-F₁₁^2-F₂₂^2-2*F₁₂


        p = λ*((F₁₁^2+F₂₁^2+F₁₂^2+F₂₂^2)-2)
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                ΔI₁ = 2*(B₁[i]*F₁₁+B₂[i]*F₂₂)
                ΔI₂ = 2*I₁*ΔI₁-2*B₁[i]*F₁₁-2*B₂[i]*F₂₂-2*B₂[i]*F₁₂
                k[2*I-1,2*J-1] +=((B₁[i]*F₁₁*B₁[j]*F₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +  (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ
                               +   B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂
                               -   p*B₁[j]*(ΔI₁*F₁₁+ΔI₁*F₁₂+(I₁-F₁₁-F₂₂-2*F₁₂)*B₂[i]+(I₁-2*F₁₁-2*F₁₂)*B₁[i]-ΔI₂))*𝑤
                              
                k[2*I-1,2*J]   += ((B₁[i]*F₁₁*B₁[j]*F₂₁ + B₂[i]*F₁₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ
                               -   p*B₂[j]*(ΔI₁*F₁₁+ΔI₁*F₁₂+(I₁-F₁₁-F₂₂-2*F₁₂)*B₂[i]+(I₁-2*F₁₁-2*F₁₂)*B₁[i]-ΔI₂))*𝑤
                               
                k[2*I,2*J-1]   += ((B₁[i]*F₂₁*B₁[j]*F₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂)*Cᵢᵢᵢᵢ
                               +   (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*Cᵢⱼᵢⱼ
                               -    p*B₁[j]*(ΔI₁*F₂₂+ΔI₁*F₁₂+(I₁-F₁₁-F₂₂-2*F₁₂)*B₁[i]+(I₁-2*F₂₂-2*F₁₂)*B₂[i]-ΔI₂))*𝑤

                k[2*I,2*J]     += ((B₁[i]*F₂₁*B₁[j]*F₂₁ + B₂[i]*F₂₂*B₂[j]*F₂₂)*Cᵢᵢᵢᵢ 
                               +   (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*Cᵢᵢⱼⱼ
                               +   (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*Cᵢⱼᵢⱼ
                               +    B₁[i]*B₁[j]*S₁₁+B₂[i]*B₂[j]*S₂₂+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂
                               -   p*B₂[j]*(ΔI₁*F₂₂+ΔI₁*F₁₂+(I₁-F₁₁-F₂₂-2*F₁₂)*B₁[i]+(I₁-2*F₂₂-2*F₁₂)*B₂[i]-ΔI₂))*𝑤

            end
        end
    end
end

function (op::Operator{:∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        Cᵢⱼᵢⱼ = E/(1+ν)/2
       
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        J = F₁₁*F₂₂-F₁₂*F₂₁
        F₁₁ = J^(-1/3)*F₁₁
        F₂₂ = J^(-1/3)*F₂₂
        F₁₂ = J^(-1/3)*F₁₂
        F₂₁ = J^(-1/3)*F₂₁

        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
function (op::Operator{:∫∫EᵢⱼSᵢⱼdxdy_SaintVenantKirchhoff_mix})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    ν = op.ν
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        # Cᵢᵢᵢᵢ = E*(1-ν)/(1+ν)/(1-2*ν)
        # Cᵢᵢⱼⱼ = E*ν/(1+ν)/(1-2*ν)
        # Cᵢⱼᵢⱼ = E/(1+ν)/2

        λ = E*ν/(1+ν)/(1-2*ν)
        μ = E/(1+ν)/2
        Cᵢᵢᵢᵢ = 2*μ
        Cᵢⱼᵢⱼ = μ
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        J = F₁₁*F₂₂-F₁₂*F₂₁
        # F₁₁ = J^(-1/3)*F₁₁
        # F₂₂ = J^(-1/3)*F₂₂
        # F₁₂ = J^(-1/3)*F₁₂
        # F₂₁ = J^(-1/3)*F₂₁

        E₁₁ = 0.5*(F₁₁*F₁₁+F₂₁*F₂₁-1.0)
        E₁₂ = 0.5*(F₁₁*F₁₂+F₂₁*F₂₂)
        E₂₂ = 0.5*(F₁₂*F₁₂+F₂₂*F₂₂-1.0)
        # S₁₁ = Cᵢᵢᵢᵢ*E₁₁+Cᵢᵢⱼⱼ*E₂₂
        # S₂₂ = Cᵢᵢⱼⱼ*E₁₁+Cᵢᵢᵢᵢ*E₂₂
        # S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        
        S₁₁ = Cᵢᵢᵢᵢ*E₁₁
        S₂₂ = Cᵢᵢᵢᵢ*E₂₂
        S₁₂ = 2.0*Cᵢⱼᵢⱼ*E₁₂
        I₁ = F₁₁*F₁₁+F₂₂*F₂₂
        I₂ = I₁^2-F₁₁^2-F₂₂^2-2*F₁₂


        p = λ*((F₁₁^2+F₂₁^2+F₁₂^2+F₂₂^2)-2)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂
                        -p*(B₁[i]+B₂[i])*(-F₁₁^2+(I₁-F₁₂)*F₁₁+I₁*F₁₂-F₁₂^2-F₁₂*F₂₂-I₂))*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂
                        -p*(B₁[i]+B₂[i])*(-F₂₂^2+(I₁-F₁₂)*F₂₂+I₁*F₁₂-F₁₂*F₁₁-F₁₂^2-I₂))*𝑤
        end
    end
end

function (op::Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
function (op::Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean_modified})(aᵤ::T,aₚ::S;k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        N = ξₚ[:𝝭]
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        J = F₁₁*F₂₂-F₁₂*F₂₁
        # F₁₁ = J^(-1/3)*F₁₁
        # F₂₂ = J^(-1/3)*F₂₂
        # F₁₂ = J^(-1/3)*F₁₂
        # F₂₁ = J^(-1/3)*F₂₁
        F₁₁ = (1/J)^(1/3)*F₁₁
        F₂₂ = (1/J)^(1/3)*F₂₂
        F₁₂ = (1/J)^(1/3)*F₁₂
        F₂₁ = (1/J)^(1/3)*F₂₁

        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0

        
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
        C⁻¹₃₃=1.0/detC*(1-I₁+I₂)
       
        
        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
        S₃₃ = λ*J*(J-1)*C⁻¹₃₃ + μ*(1-C⁻¹₃₃)
        

        # S₁₁ = μ*(1-C⁻¹₁₁)
        # S₂₂ = μ*(1-C⁻¹₂₂)
        # S₁₂ = - μ*C⁻¹₁₂
        # S₃₃ = μ*(1-C⁻¹₃₃)
        
        σ₁₁ = 1/J*F₁₁*S₁₁*F₁₁
        σ₂₂ = 1/J*F₂₂*S₂₂*F₂₂
        σ₁₂ = 1/J*F₁₂*S₁₂*F₁₂ 
        σ₃₃ = 1/J*S₃₃

        p = (σ₁₁+σ₂₂+σ₃₃)/3

        for (i,xᵢ) in  enumerate(𝓒ₚ)  
            pʰ += N[i]*xᵢ.q
    
        end
        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₃₃₃₃=λ*J*(2*J-1.0)*C⁻¹₃₃*C⁻¹₃₃+(μ- λ*J*(J-1.0))*2*C⁻¹₃₃*C⁻¹₃₃
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂) 

        C₁₁₁₁=(4*C₁₁₁₁-4*C₁₁₂₂+C₃₃₃₃+C₂₂₂₂)/9+(-4*σ₁₁+p+4*pʰ)/3
        C₂₂₂₂=(C₁₁₁₁-4*C₁₁₂₂+C₃₃₃₃+4*C₂₂₂₂)/9+(-4*σ₂₂+p+4*pʰ)/3
        C₁₁₂₂=(-2*C₁₁₁₁+5*C₁₁₂₂+C₃₃₃₃-2*C₂₂₂₂)/9+(-2*σ₁₁-2*σ₂₂+7*p+2*pʰ)/3
        C₁₁₁₂=(2*C₁₁₁₂-C₂₂₁₂-2*σ₁₂)/3
        C₂₂₁₂=(-C₁₁₁₂+2*C₂₂₁₂-2*σ₁₂)/3    
        C₁₂₁₂=(-C₁₁₁₂-C₂₂₁₂-2*σ₁₂)/3  
        for (i,xᵢ) in enumerate(𝓒ᵤ)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function (op::Operator{:Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂- λ*J*(J-1.0)*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂- λ*J*(J-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function (op::Operator{:Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        C₁₁₁₁=μ*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=μ*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=μ*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=μ*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=μ*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=μ*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function (op::Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function (op::Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean_modified})(aᵤ::T,aₚ::S;f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒
    𝓒ₚ = aₚ.𝓒
    𝓖ᵤ = aᵤ.𝓖
    𝓖ₚ = aₚ.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for (ξᵤ,ξₚ) in zip(𝓖ᵤ,𝓖ₚ)
        B₁ = ξᵤ[:∂𝝭∂x]
        B₂ = ξᵤ[:∂𝝭∂y]
        𝑤 = ξᵤ.𝑤
        N = ξₚ[:𝝭]
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒ᵤ)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end

        J = F₁₁*F₂₂-F₁₂*F₂₁
        # println(J)
        # F₁₁ = J^(-1/3)*F₁₁
        # F₂₂ = J^(-1/3)*F₂₂
        # F₁₂ = J^(-1/3)*F₁₂
        # F₂₁ = J^(-1/3)*F₂₁

        F₁₁ = (1/J)^(1/3)*F₁₁
        F₂₂ = (1/J)^(1/3)*F₂₂
        F₁₂ = (1/J)^(1/3)*F₁₂
        F₂₁ = (1/J)^(1/3)*F₂₁

        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
        C⁻¹₃₃=1.0/detC*(1-I₁+I₂)
        # S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        # S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        # S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁ + μ*(1-C⁻¹₁₁)
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂ + μ*(1-C⁻¹₂₂)
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂ - μ*C⁻¹₁₂
        S₃₃ = λ*J*(J-1)*C⁻¹₃₃ + μ*(1-C⁻¹₃₃)
        

        # S₁₁ = μ*(1-C⁻¹₁₁)
        # S₂₂ = μ*(1-C⁻¹₂₂)
        # S₁₂ = - μ*C⁻¹₁₂
        # S₃₃ = μ*(1-C⁻¹₃₃)
        
        σ₁₁ = 1/J*F₁₁*S₁₁*F₁₁
        σ₂₂ = 1/J*F₂₂*S₂₂*F₂₂
        σ₁₂ = 1/J*F₁₂*S₁₂*F₁₂ 
        σ₃₃ = 1/J*S₃₃

        p = (σ₁₁+σ₂₂+σ₃₃)/3

        for (i,xᵢ) in  enumerate(𝓒ₚ)  
            pʰ += N[i]*xᵢ.q
    
        end
        C₁₁₁₁=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₁+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₃₃₃₃=λ*J*(2*J-1.0)*C⁻¹₃₃*C⁻¹₃₃+(μ- λ*J*(J-1.0))*2*C⁻¹₃₃*C⁻¹₃₃
        C₁₁₂₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₂₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=λ*J*(2*J-1.0)*C⁻¹₁₁*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=λ*J*(2*J-1.0)*C⁻¹₂₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*2*C⁻¹₁₂*C⁻¹₂₂     
        C₁₂₁₂=λ*J*(2*J-1.0)*C⁻¹₁₂*C⁻¹₁₂+(μ- λ*J*(J-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂) 

        C₁₁₁₁=(4*C₁₁₁₁-4*C₁₁₂₂+C₃₃₃₃+C₂₂₂₂)/9+(-4*σ₁₁+p+4*pʰ)/3
        C₂₂₂₂=(C₁₁₁₁-4*C₁₁₂₂+C₃₃₃₃+4*C₂₂₂₂)/9+(-4*σ₂₂+p+4*pʰ)/3
        C₁₁₂₂=(-2*C₁₁₁₁+5*C₁₁₂₂+C₃₃₃₃-2*C₂₂₂₂)/9+(-2*σ₁₁-2*σ₂₂+7*p+2*pʰ)/3
        C₁₁₁₂=(2*C₁₁₁₂-C₂₂₁₂-2*σ₁₂)/3
        C₂₂₁₂=(-C₁₁₁₂+2*C₂₂₁₂-2*σ₁₂)/3    
        C₁₂₁₂=(-C₁₁₁₂-C₂₂₁₂-2*σ₁₂)/3  
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function (op::Operator{:∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = λ*J*(J-1)*C⁻¹₁₁
        S₂₂ = λ*J*(J-1)*C⁻¹₂₂
        S₁₂ = λ*J*(J-1)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function (op::Operator{:∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        S₁₁ = μ*(1-C⁻¹₁₁)
        S₂₂ = μ*(1-C⁻¹₂₂)
        S₁₂ = - μ*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end

function (op::Operator{:∫vᵢuᵢds})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        u₁ = 0.0
        u₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
        end
        Δu₁ = g₁-u₁
        Δu₂ = g₂-u₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += α*N[i]*(n₁₁*Δu₁+n₁₂*Δu₂)*𝑤
            f[2*I]   += α*N[i]*(n₁₂*Δu₁+n₂₂*Δu₂)*𝑤
        end
    end
end
function (op::Operator{:Δ∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        # println(J⁻²)
        J⁻²³=cbrt(J⁻²)
        # println(J⁻²³)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₁
            #   println(C₁₁₁₁)
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ+(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(J^2*K+2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂+(1/3*μ*J⁻²³*I₁-0.5*K*(J^2-1.0))*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)+0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)+0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)+0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end

function (op::Operator{:Δ∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)

        C₁₁₁₁=J^2*K*C⁻¹₁₁*C⁻¹₁₁+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=J^2*K*C⁻¹₂₂*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=J^2*K*C⁻¹₁₁*C⁻¹₂₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=J^2*K*C⁻¹₁₁*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=J^2*K*C⁻¹₂₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=J^2*K*C⁻¹₁₂*C⁻¹₁₂+0.5*K*(J^2-1.0)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
function (op::Operator{:Δ∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        C₁₁₁₁=(-2/3*C⁻¹₁₁-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₁ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₁
        C₂₂₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₂₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₂₂
        C₁₁₂₂=(-2/3*C⁻¹₂₂-2/3*C⁻¹₁₁)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₂₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₂*C⁻¹₁₂
        C₁₁₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₁₁*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₁₁*C⁻¹₁₂ 
        C₂₂₁₂=(-2/3*C⁻¹₁₂)*J⁻²³*μ + (2/9*J⁻²³*I₁*μ)*C⁻¹₂₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*2*C⁻¹₂₂*C⁻¹₁₂     
        C₁₂₁₂=(2/9*J⁻²³*I₁*μ)*C⁻¹₁₂*C⁻¹₁₂ + (1/3*μ*J⁻²³*I₁)*(C⁻¹₁₁*C⁻¹₂₂+C⁻¹₁₂*C⁻¹₁₂)   
              
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] +=(B₁[i]*F₁₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               + (B₁[i]*F₁₁*B₂[j]*F₁₂ + B₂[i]*F₁₂*B₁[j]*F₁₁)*C₁₁₂₂
                               + (B₁[i]*F₁₁*B₂[j]*F₁₁ + B₁[i]*F₁₁*B₁[j]*F₁₂)*C₁₁₁₂
                               + (B₁[i]*F₁₂*B₁[j]*F₁₁ + B₂[i]*F₁₁*B₁[j]*F₁₁)*C₁₁₁₂
                               + (B₂[i]*F₁₂*B₂[j]*F₁₁ + B₂[i]*F₁₂*B₁[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂*B₂[j]*F₁₂ + B₂[i]*F₁₁*B₂[j]*F₁₂)*C₂₂₁₂
                               + (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂
                               +  B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤
                              
                k[2*I-1,2*J]   += (B₁[i]*F₁₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₁₂*B₂[j]*F₂₂*C₂₂₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₂ + B₂[i]*F₁₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₁₁*B₂[j]*F₂₁ + B₁[i]*F₁₁*B₁[j]*F₂₂)*C₁₁₁₂ 
                               +  (B₁[i]*F₁₂*B₁[j]*F₂₁ + B₂[i]*F₁₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₁₂*B₂[j]*F₂₁ + B₂[i]*F₁₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂*B₂[j]*F₂₂ + B₂[i]*F₁₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₁₂+B₂[i]*F₁₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂)*𝑤
                               
                k[2*I,2*J-1]   += (B₁[i]*F₂₁*B₁[j]*F₁₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₁₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₂ + B₂[i]*F₂₂*B₁[j]*F₁₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₁₁ + B₁[i]*F₂₁*B₁[j]*F₁₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₁₁ + B₂[i]*F₂₁*B₁[j]*F₁₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₁₁ + B₂[i]*F₂₂*B₁[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₁₂ + B₂[i]*F₂₁*B₂[j]*F₁₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₁₂+B₂[j]*F₁₁)*C₁₂₁₂)*𝑤

                k[2*I,2*J]     += (B₁[i]*F₂₁*B₁[j]*F₂₁*C₁₁₁₁ + B₂[i]*F₂₂*B₂[j]*F₂₂*C₂₂₂₂ 
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₂ + B₂[i]*F₂₂*B₁[j]*F₂₁)*C₁₁₂₂
                               +  (B₁[i]*F₂₁*B₂[j]*F₂₁ + B₁[i]*F₂₁*B₁[j]*F₂₂)*C₁₁₁₂
                               +  (B₁[i]*F₂₂*B₁[j]*F₂₁ + B₂[i]*F₂₁*B₁[j]*F₂₁)*C₁₁₁₂
                               +  (B₂[i]*F₂₂*B₂[j]*F₂₁ + B₂[i]*F₂₂*B₁[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂*B₂[j]*F₂₂ + B₂[i]*F₂₁*B₂[j]*F₂₂)*C₂₂₁₂
                               +  (B₁[i]*F₂₂+B₂[i]*F₂₁)*(B₁[j]*F₂₂+B₂[j]*F₂₁)*C₁₂₁₂
                               +   B₁[i]*B₁[j]*S₁₁+(B₁[i]*B₂[j]+B₂[i]*B₁[j])*S₁₂+B₂[i]*B₂[j]*S₂₂)*𝑤

            end
        end
    end
end
function (op::Operator{:∫∫EᵢⱼSᵢⱼdxdy_NeoHookean2})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)
        
        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁) + 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂) + 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂) + 0.5*K*(J^2-1.0)*C⁻¹₁₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
function (op::Operator{:∫∫EᵛᵢⱼSᵛᵢⱼdxdy_NeoHookean2})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        
        S₁₁ = 0.5*K*(J^2-1.0)*C⁻¹₁₁
        S₂₂ = 0.5*K*(J^2-1.0)*C⁻¹₂₂
        S₁₂ = 0.5*K*(J^2-1.0)*C⁻¹₁₂
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end
function (op::Operator{:∫∫EᵈᵢⱼSᵈᵢⱼdxdy_NeoHookean2})(ap::T;f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E=op.E
    ν=op.ν
    K=E/(1-2*ν)/3
    λ=E*ν/(1+ν)/(1-2*ν)
    μ=E/(1+ν)/2
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        F₁₁ = 1.0
        F₁₂ = 0.0
        F₂₁ = 0.0
        F₂₂ = 1.0
        for (i,xᵢ) in  enumerate(𝓒)  
            F₁₁ += B₁[i]*xᵢ.d₁
            F₁₂ += B₂[i]*xᵢ.d₁
            F₂₁ += B₁[i]*xᵢ.d₂
            F₂₂ += B₂[i]*xᵢ.d₂
        end
        C₁₁ = F₁₁*F₁₁+F₂₁*F₂₁       
        C₁₂ = F₁₁*F₁₂+F₂₁*F₂₂
        C₂₂ = F₁₂*F₁₂+F₂₂*F₂₂
        C₃₃ = 1.0
        J = F₁₁*F₂₂-F₁₂*F₂₁
        detC = (C₁₁*C₂₂-C₁₂*C₁₂)
        I₁ = C₁₁+C₂₂+C₃₃
        I₂ = 0.5*(I₁^2-C₁₁^2-C₂₂^2-C₃₃^2-2*C₁₂^2)
        C⁻¹₁₁=1.0/detC*(C₁₁*C₁₁+C₁₂*C₁₂-I₁*C₁₁+I₂)
        C⁻¹₁₂=1.0/detC*(C₁₁*C₁₂+C₁₂*C₂₂-I₁*C₁₂)
        C⁻¹₂₂=1.0/detC*(C₁₂*C₁₂+C₂₂*C₂₂-I₁*C₂₂+I₂)

        J⁻²=(1.0/J)^2
        J⁻²³=cbrt(J⁻²)
        S₁₁ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₁₁)
        S₂₂ = μ*J⁻²³*(1.0-1/3*I₁*C⁻¹₂₂)
        S₁₂ = μ*J⁻²³*(-1/3*I₁*C⁻¹₁₂)
  
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[2*I-1] += (B₁[i]*F₁₁*S₁₁+B₂[i]*F₁₂*S₂₂+(B₁[i]*F₁₂+B₂[i]*F₁₁)*S₁₂)*𝑤
            f[2*I]   += (B₁[i]*F₂₁*S₁₁+B₂[i]*F₂₂*S₂₂+(B₁[i]*F₂₂+B₂[i]*F₂₁)*S₁₂)*𝑤
        end
    end
end