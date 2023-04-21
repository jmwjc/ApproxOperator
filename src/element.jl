"""
Element{T}<:AbstractElement{T}
"""
struct Element{T}<:AbstractElement{T}
    𝓒::Tuple{Int,Int,Vector{Node{(:𝐼,),1}}}
    𝓖::Tuple{Int,Int,Vector{Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}}}
end

# function get𝓒(a::T) where T<:AbstractElement
#     𝓒 = getfield(a,:𝓒)
#     return (𝓒[3][a.𝓒[1]+i] for i in 1:𝓒[2])
# end
# function get𝓖(a::T) where T<:AbstractElement
#     𝓖 = getfield(a,:𝓖)
#     return (𝓖[3][𝓖[1]+i] for i in 1:𝓖[2])
# end

function Base.getproperty(a::T,s::Symbol) where T<:AbstractElement
    if s∈(:𝓒,:𝓖)
        𝓐 =  getfield(a,s)
        return (𝓐[3][𝓐[1]+i] for i in 1:𝓐[2])
    else
        𝓖 = getfield(a,:𝓖)
        ξ = 𝓖[3][𝓖[1]+1]
        return getproperty(ξ,s)
    end
end

function Base.setproperty!(ap::T,s::Symbol) where T<:AbstractElement
    𝓖 = getfield(ap,:𝓖)
    ξ = 𝓖[3][𝓖[1]+1]
    setproperty!(ξ,s)
end

for set𝝭 in (:set𝝭!,:set∇𝝭!)
    @eval begin
        function $set𝝭(a::T) where T<:AbstractElement
            𝓖 = a.𝓖
            for x in 𝓖
                $set𝝭(a,x)
            end
        end
    end
end