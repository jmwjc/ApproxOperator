

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
        D¹¹¹¹ = a¹¹*a¹¹
        D²²²² = a²²*a²²
        D¹¹²² = ν*a¹¹*a²² + (1-ν)*a¹²*a¹²
        D¹¹¹² = a¹¹*a¹²
        D²²¹² = a²²*a¹²
        D¹²¹² = 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Dᵐ*(D¹¹¹¹*B₁[i]*B₁[j] + D¹¹¹²*(B₁[i]*B₂[j] + B₂[i]*B₁[j]) + D¹²¹²*B₂[i]*B₂[j])
                               +   Dᵇ*())*𝑤
                k[2*I-1,2*J]   += (Dᵐ*(D¹¹¹²*B₁[i]*B₁[j] + D¹¹²²*B₁[i]*B₂[j] + D¹²¹²*B₂[i]*B₁[j] + D²²¹²*B₂[i]*B₂[j])
                               +   Dᵇ*())*𝑤
                k[2*I,2*J-1]   += (Dᵐ*(D¹¹¹²*B₁[i]*B₁[j] + D¹²¹²*B₁[i]*B₂[j] + D¹¹²²*B₂[i]*B₁[j] + D²²¹²*B₂[i]*B₂[j])
                               +   Dᵇ*())*𝑤
                k[2*I,2*J]     += (Dᵐ*(D¹²¹²*B₁[i]*B₁[j] + D²²¹²*(B₁[i]*B₂[j] + B₂[i]*B₁[j]) + D²²²²*B₂[i]*B₂[j])
                               +   Dᵇ*())*𝑤
            end
        end
    end
end