
function set𝝭!(::Element{:Tri3},x::Node)
    𝝭 = x[:𝝭]
    𝝭[1] = x.ξ
    𝝭[2] = x.η
    𝝭[3] = 1.0-x.ξ-x.η
end
function set∇𝝭!(ap::Element{:Tri3},x::Node)
    𝐴 = ap.𝐴
    v₁,v₂,v₃ = ap.𝓒
    x₁ = v₁.x
    y₁ = v₁.y
    x₂ = v₂.x
    y₂ = v₂.y
    x₃ = v₃.x
    y₃ = v₃.y
    𝝭 = x[:𝝭]
    𝝭[1] = x.ξ
    𝝭[2] = x.η
    𝝭[3] = 1.0-x.ξ-x.η
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (y₂-y₃)/2.0/𝐴
    ∂𝝭∂x[2] = (y₃-y₁)/2.0/𝐴
    ∂𝝭∂x[3] = (y₁-y₂)/2.0/𝐴
    ∂𝝭∂y[1] = (x₃-x₂)/2.0/𝐴
    ∂𝝭∂y[2] = (x₁-x₃)/2.0/𝐴
    ∂𝝭∂y[3] = (x₂-x₁)/2.0/𝐴
end