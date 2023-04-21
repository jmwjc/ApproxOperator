
function set𝝭ₕ₁!(ap::ReproducingKernel{:Linear1D},x::Node)
    𝓒 = ap.𝓒
    𝝭 = x[:𝝭]
    𝝭₁ = x[:𝝭₁]
    𝗠⁻¹ = cal𝗠!(ap,x)
    for (i,xᵢ) in enumerate(𝓒)
        r = ((x.x-xᵢ.x)^2+(x.y-xᵢ.y)^2+(x.z-xᵢ.z)^2)^0.5
        𝝭_ = 0.0
        𝝭₁_ = 0.0
        if r < xᵢ.s
            m₀ = xᵢ.m₀
            m₁ = xᵢ.m₁
            𝒑ᵢ,𝒑ᵢ₁ = get∇₁𝒑(ap,xᵢ)
            𝒑 = get𝒑(ap,x)
            𝗠⁻¹𝒑 = 𝗠⁻¹*𝒑
            for j in 1:2
                𝝭_  += (m₀*𝒑ᵢ[j]+m₁*𝒑ᵢ₁[j])*𝗠⁻¹𝒑[j]
                𝝭₁_ += m₀*𝒑ᵢ₁[j]*𝗠⁻¹𝒑[j]
            end
        end
        𝝭[i] = 𝝭_
        𝝭₁[i] = 𝝭₁_
    end
end