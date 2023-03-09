
@inline get𝒙(ap::T,ξ::Node) where T<:AbstractElement{:Tri6} = get𝒙(ap,ξ.ξ,ξ.η)

function get𝒙(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Tri6}
    γ = 1.0-ξ-η
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y;z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y;z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y;z₃ = ap.𝓒[3].z
    x₄ = ap.𝓒[4].x;y₄ = ap.𝓒[4].y;z₄ = ap.𝓒[4].z
    x₅ = ap.𝓒[5].x;y₅ = ap.𝓒[5].y;z₅ = ap.𝓒[5].z
    x₆ = ap.𝓒[6].x;y₆ = ap.𝓒[6].y;z₆ = ap.𝓒[6].z
    N₁ = ξ*(2*ξ-1)
    N₂ = η*(2*η-1)
    N₃ = γ*(2*γ-1)
    N₄ = 4*ξ*η
    N₅ = 4*η*γ
    N₆ = 4*γ*ξ
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄+x₅*N₅+x₆*N₆,
            y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄+y₅*N₅+y₆*N₆,
            z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄+z₅*N₅+z₆*N₆)
end

@inline get𝑤(ap::T,ξ::Node) where T<:AbstractElement{:Tri6} = get𝐴(ap)*ξ.w

function get𝐴(ap::T) where T<:AbstractElement{:Tri6}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    𝐴₁ = 0.5*(y₁*z₂+y₂*z₃+y₃*z₁-y₂*z₁-y₃*z₂-y₁*z₃)
    𝐴₂ = 0.5*(z₁*x₂+z₂*x₃+z₃*x₁-z₂*x₁-z₃*x₂-z₁*x₃)
    𝐴₃ = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
    return (𝐴₁^2 + 𝐴₂^2 + 𝐴₃^2)^0.5
end

function set𝝭!(ap::Element{:Tri6},x::Node)
    𝝭 = x[:𝝭]
    ξ = x.ξ
    η = x.η
    γ = 1.0-ξ-η
    𝝭[1] = ξ*(2*ξ-1)
    𝝭[2] = η*(2*η-1)
    𝝭[3] = γ*(2*γ-1)
    𝝭[4] = 4*ξ*η
    𝝭[5] = 4*η*γ
    𝝭[6] = 4*γ*ξ
end

function set∇𝝭!(ap::Element{:Tri6},x::Node)
    𝐴 = get𝐴(ap)
    ξ = x.ξ
    η = x.η
    γ = 1.0-ξ-η
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    J₁₁ = (y₂-y₃)/2.0/𝐴
    J₂₁ = (y₃-y₁)/2.0/𝐴
    J₃₁ = (y₁-y₂)/2.0/𝐴
    J₁₂ = (x₃-x₂)/2.0/𝐴
    J₂₂ = (x₁-x₃)/2.0/𝐴
    J₃₂ = (x₂-x₁)/2.0/𝐴
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