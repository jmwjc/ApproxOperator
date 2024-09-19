"""
Operator
"""
struct Operator{T}
    data::Dict{Symbol,Float64}
end
Operator{T}(d::Pair{Symbol}...) where T = Operator{T}(Dict(d))


getproperty(op::Operator,f::Symbol) = getfield(op,:data)[f]

(op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},g::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,k=k,g=g)
(op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,k=k,f=f)
# (op::Operator)(aps::Vector{T},gps::Vector{S},g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,g=g,q=q)
(op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,k=k)
(op::Operator)(aps::Vector{T},gps::Vector{S},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement} = op.(aps,gps,f=f)
(op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement = op.(aps,k=k,f=f)
(op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64}) where T<:AbstractElement = op.(aps,k=k)
(op::Operator)(aps::Vector{T},f::AbstractVector{Float64}) where T<:AbstractElement = op.(aps,f=f)
(op::Operator)(aps::Vector{T}) where T<:AbstractElement = op.(aps)

function prescribe!(ap::T,sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    𝓖 = ap.𝓖
    s,f = sf
    for ξ in 𝓖
        𝒙 = (ξ.x,ξ.y,ξ.z)
        if applicable(f,𝒙...)
            v = f(𝒙...)
        elseif applicable(f,𝒙...,ξ.n₁)
            v = f(𝒙...,ξ.n₁)
        elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂)
            v = f(𝒙...,ξ.n₁,ξ.n₂)
        elseif applicable(f,𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
            v = f(𝒙...,ξ.n₁,ξ.n₂,ξ.n₃)
        end
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},sf::Pair{Symbol,F}) where {T<:AbstractElement,F<:Function}
    s,f = sf
    for ap in aps
        data = getfield((ap.𝓖)[1],:data)
        n = length(data[:x][2])
        haskey(data,s) ? nothing : push!(data,s=>(2,zeros(n)))
        prescribe!(ap,sf)
    end
end