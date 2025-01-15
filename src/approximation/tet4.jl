
function set𝝭!(::Element{:Tet4},x::Node)
    ξ = x.ξ
    η = x.η
    γ = x.γ
    𝝭 = x[:𝝭]
   
    𝝭[1] = 1.0-ξ-η-γ
    𝝭[2] = ξ
    𝝭[3] = η
    𝝭[4] = γ
end
function set∇𝝭!(ap::Element{:Tet4},x::Node)
    𝐽 = ap.𝐽
    ξ = x.ξ
    η = x.η
    γ = x.γ
    𝝭 = x[:𝝭]
    𝝭[1] = 1.0-ξ-η-γ
    𝝭[2] = ξ
    𝝭[3] = η
    𝝭[4] = γ
    𝓒 = ap.𝓒
    v₁,v₂,v₃,v₄ = ap.𝓒
    x₁ = v₁.x
    y₁ = v₁.y
    z₁ = v₁.z
    x₂ = v₂.x
    y₂ = v₂.y
    z₂ = v₂.z
    x₃ = v₃.x
    y₃ = v₃.y
    z₃ = v₃.z
    x₄ = v₄.x
    y₄ = v₄.y
    z₄ = v₄.z
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂z = x[:∂𝝭∂z]
    ∂𝝭∂x[1] = (y₃*z₂+y₄*z₃+y₂*z₄-y₂*z₃-y₃*z₄-y₄*z₂)/𝐽
    ∂𝝭∂x[2] = (y₄*z₁+y₃*z₄+y₁*z₃-y₁*z₄-y₄*z₃-y₃*z₁)/𝐽
    ∂𝝭∂x[3] = (y₂*z₁+y₄*z₂+y₁*z₄-y₁*z₂-y₂*z₄-y₄*z₁)/𝐽
    ∂𝝭∂x[4] = (y₃*z₁+y₂*z₃+y₁*z₂-y₁*z₃-y₃*z₂-y₂*z₁)/𝐽
    ∂𝝭∂y[1] = (z₃*x₂+z₄*x₃+z₂*x₄-z₂*x₃-z₃*x₄-z₄*x₂)/𝐽
    ∂𝝭∂y[2] = (z₄*x₁+z₃*x₄+z₁*x₃-z₁*x₄-z₄*x₃-z₃*x₁)/𝐽
    ∂𝝭∂y[3] = (z₂*x₁+z₄*x₂+z₁*x₄-z₁*x₂-z₂*x₄-z₄*x₁)/𝐽
    ∂𝝭∂y[4] = (z₃*x₁+z₂*x₃+z₁*x₂-z₁*x₃-z₃*x₂-z₂*x₁)/𝐽
    ∂𝝭∂z[1] = (x₃*y₂+x₄*y₃+x₂*y₄-x₂*y₃-x₃*y₄-x₄*y₂)/𝐽
    ∂𝝭∂z[2] = (x₄*y₁+x₃*y₄+x₁*y₃-x₁*y₄-x₄*y₃-x₃*y₁)/𝐽
    ∂𝝭∂z[3] = (x₂*y₁+x₄*y₂+x₁*y₄-x₁*y₂-x₂*y₄-x₄*y₁)/𝐽
    ∂𝝭∂z[4] = (x₃*y₁+x₂*y₃+x₁*y₂-x₁*y₃-x₃*y₂-x₂*y₁)/𝐽
end