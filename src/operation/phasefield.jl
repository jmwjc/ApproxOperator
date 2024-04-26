
function (op::Operator{:∫v²uₓuₓdx})(ap::T;k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*B[i]*B[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx_hard_device})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
        end
        ℋₜ = max(ℋ,(ε-1)^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*kc/2/l*𝑤
        end
    end
end

function (op::Operator{:∫vₓvₓvvdx})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        ε = 0.0
        v = 0.0
        ∂v∂x = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
            ∂v∂x += B[i]*xᵢ.v
        end
        ℋₜ = max(ℋ,(v^2+η)*ε^2)
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*B[i]*B[j] + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += N[i]*(kc/2/l - η*ℋₜ)*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_1D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε += B[i]*xᵢ.u
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,(v^2+η)*ε^2)
    end
end
function (op::Operator{:∫∫∇v∇vvvdxdy})(ap::T;k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kc = op.k
    l = op.l
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        ℋ = ξ.ℋ
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        v = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            v += N[i]*xᵢ.v
            ∂v∂x += B₁[i]*xᵢ.v
            ∂v∂y += B₂[i]*xᵢ.v
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
        end
        εᵢᵢ = ε₁₁ + ε₂₂
        e₁₁ = ε₁₁ - εᵢᵢ/3
        e₂₂ = ε₂₂ - εᵢᵢ/3
        e₁₂ = ε₁₂
        sᵢᵢ = σ₁₁ + σ₂₂
        p = sᵢᵢ/3
        s₁₁ = σ₁₁ - p
        s₂₂ = σ₂₂ - p
        s₁₂ = σ₁₂
        #ℋₜ = max(ℋ,(e₁₁*s₁₁ + e₂₂*s₂₂ + e₁₂*s₁₂))
        ℋₜ = max(ℋ,(ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂))
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (kc*(2*l*(B₁[i]*B₁[j] + B₂[i]*B₂[j]) + N[i]*N[j]/2/l) + ℋₜ*N[i]*N[j])*𝑤
            end
            f[I] += (N[i]*kc/2/l - kc*(2*l*(B₁[i]*∂v∂x + B₂[i]*∂v∂y) + N[i]*v/2/l) - ℋₜ*N[i]*v)*𝑤
        end
    end
end

function (op::Operator{:UPDATE_PFM_2D})(ap::T) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂x]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        v = 0.0
        ε₁₁ = 0.0
        ε₂₂ = 0.0
        ε₁₂ = 0.0
        ℋ = ξ.ℋ
        for (i,xᵢ) in enumerate(𝓒)
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            v += N[i]*xᵢ.v
        end
        ξ.ℋ = max(ℋ,ε₁₁*σ₁₁ + ε₂₂*σ₂₂ + ε₁₂*σ₁₂)
    end
end

function (op::Operator{:CRACK_NORMAL})(aps::Vector{T},nodes::Vector{N},v::Vector{Float64}) where {T<:AbstractElement,N<:Node}
    l = op.l
    Γtmp = findall(x->x<0.02,v)
    Γfinal = Int[]
    # 循环，当Γtmp不为空时
    if ~isempty(Γtmp)
        while ~isempty(Γtmp)
            # 找到identity中v[Γtmp]的最小值，并将其索引N1记录
            _,N1 = findmin(identity,v[Γtmp])
            # 将Γtmp[N1]添加到Γfinal中
            push!(Γfinal,Γtmp[N1])
            # 记录Γtmp[N1]的x坐标和y坐标
            x₁ = nodes[Γtmp[N1]].x
            y₁ = nodes[Γtmp[N1]].y
            # 遍历Γtmp中的每一个元素
            for index in Γtmp
                # 记录当前元素对应的节点
                node = nodes[index]
                # 记录当前节点对应的x坐标和y坐标
                x₂ = node.x
                y₂ = node.y
                # 计算当前节点和Γtmp[N1]之间的距离
                Δ = ((x₁-x₂)^2 + (y₁-y₂)^2)^0.5
                # 如果距离小于l，则将当前节点从Γtmp中移除
                Δ < l ? setdiff!(Γtmp,index) : nothing
            end
        end
        sort!(Γfinal, by = i->nodes[i].x)
        n₁ = zeros(length(Γfinal)-1)
        n₂ = zeros(length(Γfinal)-1)
        for i in 1:length(Γfinal)-1
            i₁ = Γfinal[i]
            i₂ = Γfinal[i+1]
            x₁ = nodes[i₁].x
            y₁ = nodes[i₁].y
            x₂ = nodes[i₂].x
            y₂ = nodes[i₂].y
            𝐿 = ((x₁-x₂)^2 + (y₁-y₂)^2)^0.5
            n₁[i] = (y₂-y₁)/𝐿
            n₂[i] = (x₁-x₂)/𝐿
        end
        for ap in aps
            𝓖 = ap.𝓖
            for ξ in 𝓖
                𝐿 = Inf
                x₀ = ξ.x
                y₀ = ξ.y
                for i in 1:length(Γfinal)-1
                    i₁ = Γfinal[i]
                    i₂ = Γfinal[i+1]
                    x₁ = nodes[i₁].x
                    y₁ = nodes[i₁].y
                    x₂ = nodes[i₂].x
                    y₂ = nodes[i₂].y
                    A = y₁-y₂
                    B = x₂-x₁
                    C = y₁*x₂ - x₁*y₂
                    𝐿tmp = abs(A*x₀+B*y₀+C)/(A^2+B^2)^0.5
                    if 𝐿tmp < 𝐿
                        ξ.n₁ = n₁[i]
                        ξ.n₂ = n₂[i]
                        𝐿 = 𝐿tmp
                    end
                end
            end
        end
    end
end
