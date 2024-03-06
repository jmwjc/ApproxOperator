"""
Returen value (RV)

This is a simple struct
"""
struct RV
    i::Int
    v::Vector{Float64}
end

getindex(r::RV,i::Int) = r.v[r.i+i]
function setindex!(r::RV,x::Float64,i::Int)
    r.v[r.i+i] = x
end

"""
Node
"""
struct Node{T,N}
    index::NTuple{N,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const 𝑿ᵢ = Node{(:𝐼,),1}
const 𝑿ₛ = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}

function Base.getproperty(p::Node{T,N},s::Symbol) where {T,N}
    index = getfield(p,:index)
    if s ∈ T
        return index[findfirst(x->x==s,T)]
    else
        i,v = getfield(p,:data)[s]
        return v[index[i]]
    end
end
function Base.setproperty!(p::Node,s::Symbol,x::Float64)
    i,v = getfield(p,:data)[s]
    j = getfield(p,:index)[i]
    v[j] = x
end
function Base.getindex(p::Node,f::Symbol)
    i,v = getfield(p,:data)[f]
    j = getfield(p,:index)[i]
    return RV(j,v)
end

+(a::T,b::S) where {T<:Node,S<:Node} = (a.x+b.x,a.y+b.y,a.z+b.z)
-(a::T,b::S) where {T<:Node,S<:Node} = (a.x-b.x,a.y-b.y,a.z-b.z)

push!(ps::Vector{T},svs::Pair{Symbol,Vector{Float64}}...) where T<:Node = push!(ps[1],svs...)
function push!(p::Node,svs::Pair{Symbol,Vector{Float64}}...)
    for (s,v) in svs
        push!(getfield(p,:data),s=>(1,v))
    end
end

function printinfo(p::Node{S,N}) where {S,N}
    index = getfield(p,:index)
    data = getfield(p,:data)
    print("Node( ")
    for (n,s) in enumerate(S)
        i = index[n]
        print(string(s)*" = $i ")
    end
    shapes = Symbol[]
    println("):")
    for (name,(n,vs)) in data
        if n ≠ 0
            s = S[n]
            if s ≠ :𝑠
                i = index[n]
                v = vs[i]
                @printf "  %s(%s = %i): %e\n" string(name) string(s) i v
            else
                push!(shapes,name)
            end
        end
    end
    return shapes
end

Base.show(io::IO,::MIME"text/plain",p::Node) = printinfo(p)