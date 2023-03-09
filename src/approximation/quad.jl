
@inline get𝒙(ap::T,ξ::Node) where T<:AbstractElement{:Quad} = get𝒙(ap,ξ.ξ,ξ.η)

function get𝒙(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Quad}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    x₄ = ap.𝓒[4].x
    y₄ = ap.𝓒[4].y
    z₄ = ap.𝓒[4].z
    N₁,N₂,N₃,N₄ = get𝝭(ap,ξ,η)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄,y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄,z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄)
end

function get𝑱(ap::T,ξ::Node) where T<:AbstractElement{:Quad}
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(ap,ξ)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(ap,ξ)
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    return J₁₁,J₂₁,J₁₂,J₂₂
end

@inline function get𝐽(ap::T,ξ::Node) where T<:AbstractElement{:Quad}
    J₁₁,J₂₁,J₁₂,J₂₂ = get𝑱(ap,ξ)
    return J₁₁*J₂₂-J₂₁*J₁₂
end

@inline get𝑤(ap::T,ξ::Node) where T<:AbstractElement{:Quad} = get𝐽(ap,ξ)*ξ.w

function set𝝭!(ap::Element{:Quad},x::Node)
    ξ = x.ξ
    η = x.η
    𝝭 = x[:𝝭]
    𝝭[1] = 0.25*(1.0-ξ)*(1.0-η)
    𝝭[2] = 0.25*(1.0+ξ)*(1.0-η)
    𝝭[3] = 0.25*(1.0+ξ)*(1.0+η)
    𝝭[4] = 0.25*(1.0-ξ)*(1.0+η)
end

function set∇𝝭!(ap::Element{:Quad},x::Node)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(ap,x)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(ap,x)
    ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
    ∂ξ∂x =   ∂y∂η/detJ
    ∂η∂x = - ∂y∂ξ/detJ
    ∂ξ∂y = - ∂x∂η/detJ
    ∂η∂y =   ∂x∂ξ/detJ
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = ∂N₁∂ξ*∂ξ∂x + ∂N₁∂η*∂η∂x
    ∂𝝭∂x[2] = ∂N₂∂ξ*∂ξ∂x + ∂N₂∂η*∂η∂x
    ∂𝝭∂x[3] = ∂N₃∂ξ*∂ξ∂x + ∂N₃∂η*∂η∂x
    ∂𝝭∂x[4] = ∂N₄∂ξ*∂ξ∂x + ∂N₄∂η*∂η∂x
    ∂𝝭∂y[1] = ∂N₁∂ξ*∂ξ∂y + ∂N₁∂η*∂η∂y
    ∂𝝭∂y[2] = ∂N₂∂ξ*∂ξ∂y + ∂N₂∂η*∂η∂y
    ∂𝝭∂y[3] = ∂N₃∂ξ*∂ξ∂y + ∂N₃∂η*∂η∂y
    ∂𝝭∂y[4] = ∂N₄∂ξ*∂ξ∂y + ∂N₄∂η*∂η∂y
end

function get∂𝝭∂ξ(ap::Element{:Quad},η::Float64)
    ∂N₁∂ξ = - 0.25*(1-η)
    ∂N₂∂ξ =   0.25*(1-η)
    ∂N₃∂ξ =   0.25*(1+η)
    ∂N₄∂ξ = - 0.25*(1+η)
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get∂𝝭∂η(ap::Element{:Quad},ξ::Float64)
    ∂N₁∂η = - 0.25*(1-ξ)
    ∂N₂∂η = - 0.25*(1+ξ)
    ∂N₃∂η =   0.25*(1+ξ)
    ∂N₄∂η =   0.25*(1-ξ)
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end