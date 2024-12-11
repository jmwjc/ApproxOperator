# FIXME: The shape function order is wrong!
function set𝝭!(::Element{:Tri6},x::Node)
    𝝭 = x[:𝝭]
    ξ = 1.0-x.ξ-x.η
    η = x.ξ
    γ = x.η
    𝝭[1] = ξ*(2*ξ-1)
    𝝭[2] = η*(2*η-1)
    𝝭[3] = γ*(2*γ-1)
    𝝭[4] = 4*ξ*η
    𝝭[5] = 4*η*γ
    𝝭[6] = 4*γ*ξ
end

function set∇𝝭!(ap::Element{:Tri6},x::Node)
    𝐽 = x.𝐽
    ξ = 1.0-x.ξ-x.η
    η = x.ξ
    γ = x.η
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
<<<<<<< HEAD
    
    

    J₁₁ = (y₃-y₁)/𝐽
    J₂₁ = (y₁-y₂)/𝐽
    J₃₁ = (y₂-y₃)/𝐽
    J₁₂ = (x₁-x₃)/𝐽
    J₂₂ = (x₂-x₁)/𝐽
    J₃₂ = (x₃-x₂)/𝐽
=======
    J₁₁ = (y₂-y₃)/𝐽
    J₂₁ = (y₃-y₁)/𝐽
    J₃₁ = (y₁-y₂)/𝐽
    J₁₂ = (x₃-x₂)/𝐽
    J₂₂ = (x₁-x₃)/𝐽
    J₃₂ = (x₂-x₁)/𝐽
>>>>>>> f3047fcf3099bd6c9ccdda3305811c92e7a79972
    𝝭 = x[:𝝭]
    𝝭[1] = ξ*(2*ξ-1)
    𝝭[2] = η*(2*η-1)
    𝝭[3] = γ*(2*γ-1)
    𝝭[4] = 4*ξ*η
    𝝭[5] = 4*η*γ
    𝝭[6] = 4*γ*ξ
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (4.0*ξ-1.0)*J₁₁
    ∂𝝭∂x[2] = (4.0*η-1.0)*J₂₁
    ∂𝝭∂x[3] = (4.0*γ-1.0)*J₃₁
    ∂𝝭∂x[4] = 4.0*(η*J₁₁+ξ*J₂₁)
    ∂𝝭∂x[5] = 4.0*(γ*J₂₁+η*J₃₁)
    ∂𝝭∂x[6] = 4.0*(ξ*J₃₁+γ*J₁₁)
    ∂𝝭∂y[1] = (4.0*ξ-1.0)*J₁₂
    ∂𝝭∂y[2] = (4.0*η-1.0)*J₂₂
    ∂𝝭∂y[3] = (4.0*γ-1.0)*J₃₂
    ∂𝝭∂y[4] = 4.0*(η*J₁₂+ξ*J₂₂)
    ∂𝝭∂y[5] = 4.0*(γ*J₂₂+η*J₃₂)
    ∂𝝭∂y[6] = 4.0*(ξ*J₃₂+γ*J₁₂)
end