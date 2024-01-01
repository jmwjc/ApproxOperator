
function set𝝭!(::Element{:Seg2},x::Node)
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*(1.0-x.ξ)
    𝝭[2] = 0.5*(1.0+x.ξ)
end

function set∇𝝭!(ap::Element{:Seg2},x::Node)
    𝐿 = ap.𝐿
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*(1.0-x.ξ)
    𝝭[2] = 0.5*(1.0+x.ξ)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂x[1] = -1.0/𝐿
    ∂𝝭∂x[2] = 1.0/𝐿
end