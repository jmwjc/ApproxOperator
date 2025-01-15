"""
Operator
"""

function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k,f)
    end
end

function (op::Pair{F,Vector{T}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,k,f)
    end
    return k,f
end

function (ops::Vector{Pair{F,Vector{T}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(k)
    end
end

function (op::Pair{F,Vector{T}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,k)
    end
    return k
end

function (ops::Vector{Pair{F,Vector{T}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement}
    for op in ops
        op(f)
    end
end

function (op::Pair{F,Vector{T}})(f::AbstractVector) where {F<:Function,T<:AbstractElement}
    form, elms = op
    for elm in elms
        form(elm,f)
    end
    return f
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k,f)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix,f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,k,f)
    end
    return k,f
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(k)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(k::AbstractMatrix) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,k)
    end
    return k
end

function (ops::Vector{Pair{F,Tuple{Vector{T},Vector{S}}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    for op in ops
        op(f)
    end
end

function (op::Pair{F,Tuple{Vector{T},Vector{S}}})(f::AbstractVector) where {F<:Function,T<:AbstractElement,S<:AbstractElement}
    form, elms = op
    for (a,b) in zip(elms...)
        form(a,b,f)
    end
    return f
end

function getDOFs(aps::Vector{T}) where T<:AbstractElement
    𝐼 = Set{Int}()
    for ap in aps
        for xᵢ in ap.𝓒
            push!(𝐼,xᵢ.𝐼)
        end
    end
    return 𝐼
end

