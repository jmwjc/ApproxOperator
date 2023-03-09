
struct SymMat
    n::Int
    m::Vector{Float64}
end
SymMat(n::Int) = SymMat(n,zeros(Int(n*(n+1)/2)))

@inline function getindex(A::SymMat,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] : A.m[Int(i+j*(j-1)/2)]
end

@inline function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] = val : A.m[Int(i+j*(j-1)/2)] = val
end
@inline function setindex!(A::SymMat,val::Float64,i::Int)
    A.m[i] = val
end

@inline *(A::SymMat,v::NTuple{N,Float64}) where N = sum(A[1,i]*v[i] for i in 1:N)
@inline function *(v::NTuple{N,Float64},A::SymMat) where N
    for j in 1:N
        A[1,j] = sum(v[i]*A[i,j] for i in 1:N)
    end
    return A
end
@inline function -(A::SymMat)
    A.m .= .-A.m
    return A
end
@inline fill!(A::SymMat,val::Float64) = fill!(A.m,val)
function inverse!(A::SymMat)
    n = A.n
    for i in 1:n
        A[i,i] = 1.0/A[i,i]
        for j in i+1:n
            A[i,j] = - sum(A[i,k]*A[k,j] for k in i:j-1)/A[j,j]
        end
    end
    return A
end

function UUᵀ!(A::SymMat)
    n = A.n
    for i in 1:n
        for j in 1:i
            A[i,j] = sum(A[i,k]*A[k,j] for k in i:n)
        end
    end
    return A
end

function UᵀAU!(A::SymMat,U::SymMat)
    n = A.n
    for i in n:-1:1
        for j in n:-1:i
            A[i,j] = sum(U[k,i]*A[k,l]*U[l,j] for k in 1:i for l in 1:j)
        end
    end
    return A
end

function UAUᵀ!(A::SymMat,U::SymMat)
    n = A.n
    for i in 1:n
        for j in i:n
            A[i,j] = sum(U[i,k]*A[k,l]*U[j,l] for k in i:n for l in j:n)
        end
    end
    return A
end

function UUᵀAUUᵀ!(A::SymMat,U::SymMat)
    UᵀAU!(A,U)
    UAUᵀ!(A,U)
    return A
end

function cholesky!(A::SymMat)
    n = A.n
    for i in 1:n
        for k in 1:i-1
            A[i,i] -= A[k,i]^2
        end
        A[i,i] = A[i,i]^0.5
        for j in i+1:n
            for k in 1:i-1
                A[i,j] -= A[k,i]A[k,j]
            end
            A[i,j] = A[i,j]/A[i,i]
        end
    end
    return A
end

function get𝗠(ap::ReproducingKernel,s::Symbol)
    n = get𝑛𝒑(ap)
    data = getfield(ap.𝓖[1],:data)
    fill!(data[s][2],0.)
    return SymMat(n,data[s][2])
end
function get𝗚(ap::ReproducingKernel,s::Symbol)
    n = get𝑛𝒑₁(ap)
    data = getfield(ap.𝓖[1],:data)
    fill!(data[s][2],0.)
    return SymMat(n,data[s][2])
end

function cal𝗠ₕ₁₀!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        r = ((x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2)^0.5
        if r < xᵢ.s
            m₀ = xᵢ.m₀
            𝒑 = get𝒑(ap,xᵢ)
            for I in 1:n
                for J in 1:I
                    𝗠[I,J] += m₀*𝒑[I]*𝒑[J]
                end
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
    return 𝗠
end

function cal𝗠ₕ₁₁!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        r = ((x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2)^0.5
        if r < xᵢ.s
            m₀ = xᵢ.m₀
            m₁ = xᵢ.m₁
            𝒑,∂𝒑∂x = get∇₁𝒑(ap,xᵢ)
            for I in 1:n
                for J in 1:I
                    𝗠[I,J] += m₀*𝒑[I]*𝒑[J]
                            + m₁*(∂𝒑∂x[I]*𝒑[J] + 𝒑[I]*∂𝒑∂x[J])
                end
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
    return 𝗠
end

function cal𝗠ₕ₂!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        r = ((x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2 + (x.x - xᵢ.x)^2)^0.5
        if r < xᵢ.s
            m₀ = xᵢ.m₀
            m₁ = xᵢ.m₁
            m₂ = xᵢ.m₂
            𝒑,∂𝒑∂x,∂𝒑∂y = get∇₂𝒑(ap,xᵢ)
            for I in 1:n
                for J in 1:I
                    𝗠[I,J] += m₀*𝒑[I]*𝒑[J]
                            + m₁*(∂𝒑∂x[I]*𝒑[J] + 𝒑[I]*∂𝒑∂x[J])
                            + m₂*(∂𝒑∂y[I]*𝒑[J] + 𝒑[I]*∂𝒑∂y[J])
                end
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
    return 𝗠
end

function cal𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
    return 𝗠
end

function cal∇₁𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x = get∇𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x
end

function cal∇₂𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y
end

function cal∇𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
    ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠∂z,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂𝗠⁻¹∂z
end

function cal∇²₂𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]
                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]
                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end

    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y²
end

function cal∇²𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
    ∂²𝗠∂z² = get𝗠(ap,:∂²𝗠∂z²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
    ∂²𝗠∂x∂z = get𝗠(ap,:∂²𝗠∂x∂z)
    ∂²𝗠∂y∂z = get𝗠(ap,:∂²𝗠∂y∂z)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]

                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                ∂²𝗠∂z²[I,J] += ∂²𝜙∂z²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂z²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂z²[J] + 2.0*∂𝜙∂z*∂𝒑∂z[I]*𝒑[J] + 2.0*∂𝜙∂z*𝒑[I]*∂𝒑∂z[J] + 2.0*𝜙*∂𝒑∂z[I]*∂𝒑∂z[J]

                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂²𝗠∂x∂z[I,J] += ∂²𝜙∂x∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂z[J] + ∂𝜙∂x*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂x[J]

                ∂²𝗠∂y∂z[I,J] += ∂²𝜙∂y∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y∂z[J] + ∂𝜙∂y*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂𝗠∂zU = UᵀAU!(∂𝗠∂z,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂z²U = UᵀAU!(∂²𝗠∂z²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    Uᵀ∂²𝗠∂x∂zU = UᵀAU!(∂²𝗠∂x∂z,U)
    Uᵀ∂²𝗠∂y∂zU = UᵀAU!(∂²𝗠∂y∂z,U)
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂z²U[i,j] -= 2*Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂zU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂x∂zU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y∂zU[i,j] -= Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂yU[k,j]
            end
        end
    end

    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂z² = - UAUᵀ!(Uᵀ∂²𝗠∂z²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂²𝗠⁻¹∂x∂z = - UAUᵀ!(Uᵀ∂²𝗠∂x∂zU,U)
    ∂²𝗠⁻¹∂y∂z = - UAUᵀ!(Uᵀ∂²𝗠∂y∂zU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    ∂𝗠⁻¹∂z = - UAUᵀ!(Uᵀ∂𝗠∂zU,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂𝗠⁻¹∂z, ∂²𝗠⁻¹∂x∂z, ∂²𝗠⁻¹∂y∂z, ∂²𝗠⁻¹∂z²
end

function cal∇³𝗠!(ap::ReproducingKernel,x::Node)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝗠∂x²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝗠∂x∂y)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝗠∂y²)
    ∂³𝗠∂x³ = get𝗠(ap,:∂³𝗠∂x³)
    ∂³𝗠∂x²∂y = get𝗠(ap,:∂³𝗠∂x²∂y)
    ∂³𝗠∂x∂y² = get𝗠(ap,:∂³𝗠∂x∂y²)
    ∂³𝗠∂y³ = get𝗠(ap,:∂³𝗠∂y³)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]

                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                ∂³𝗠∂x³[I,J] += ∂³𝜙∂x³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x³[J] + 3.0*∂²𝜙∂x²*∂𝒑∂x[I]*𝒑[J] + 3.0*∂𝜙∂x*∂²𝒑∂x²[I]*𝒑[J] + 3.0*∂²𝜙∂x²*𝒑[I]*∂𝒑∂x[J] + 3.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x²[J] + 3.0*𝜙*∂²𝒑∂x²[I]*∂𝒑∂x[J] + 3.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x²[J] + 6.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂³𝗠∂x²∂y[I,J] += ∂³𝜙∂x²∂y*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x²∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x²∂y[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂x[I]*𝒑[J] + ∂²𝜙∂x²*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂x*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂y*∂²𝒑∂x²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂x[J] + ∂²𝜙∂x²*𝒑[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂y*𝒑[I]*∂²𝒑∂x²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂x[J] + 𝜙*∂²𝒑∂x²[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂y[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂y[J]

                ∂³𝗠∂x∂y²[I,J] += ∂³𝜙∂x∂y²*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x∂y²[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂y[I]*𝒑[J] + ∂²𝜙∂y²*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂y*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂x*∂²𝒑∂y²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂y[J] + ∂²𝜙∂y²*𝒑[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*𝒑[I]*∂²𝒑∂y²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂y[J] + 𝜙*∂²𝒑∂y²[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂x[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂³𝗠∂y³[I,J] += ∂³𝜙∂y³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂y³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂y³[J] + 3.0*∂²𝜙∂y²*∂𝒑∂y[I]*𝒑[J] + 3.0*∂𝜙∂y*∂²𝒑∂y²[I]*𝒑[J] + 3.0*∂²𝜙∂y²*𝒑[I]*∂𝒑∂y[J] + 3.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂y²[J] + 3.0*𝜙*∂²𝒑∂y²[I]*∂𝒑∂y[J] + 3.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂y²[J] + 6.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    Uᵀ∂³𝗠∂x³U = UᵀAU!(∂³𝗠∂x³,U)
    Uᵀ∂³𝗠∂x²∂yU = UᵀAU!(∂³𝗠∂x²∂y,U)
    Uᵀ∂³𝗠∂x∂y²U = UᵀAU!(∂³𝗠∂x∂y²,U)
    Uᵀ∂³𝗠∂y³U = UᵀAU!(∂³𝗠∂y³,U)

    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]+Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]+Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
            end
        end
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
                Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
            end
        end
    end

    ∂³𝗠⁻¹∂x³ = - UAUᵀ!(Uᵀ∂³𝗠∂x³U,U)
    ∂³𝗠⁻¹∂x²∂y = - UAUᵀ!(Uᵀ∂³𝗠∂x²∂yU,U)
    ∂³𝗠⁻¹∂x∂y² = - UAUᵀ!(Uᵀ∂³𝗠∂x∂y²U,U)
    ∂³𝗠⁻¹∂y³ = - UAUᵀ!(Uᵀ∂³𝗠∂y³U,U)
    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    𝗠⁻¹ = UUᵀ!(U)

    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂³𝗠⁻¹∂x³, ∂³𝗠⁻¹∂x²∂y, ∂³𝗠⁻¹∂x∂y², ∂³𝗠⁻¹∂y³
end

function cal𝗚!(ap::ReproducingKernel{𝑝,𝑠,𝜙,:Poi1}) where {𝑝,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    𝗚 = get𝗚(ap,:∇̃)
    x = 𝓒[1]
    n = get𝑛𝒑₁(ap)
    for ξ in 𝓖
        Δx = x - ξ
        𝒒 = get𝒑₁(ap,ξ)
        w = get𝜙(ap,x,Δx)
        for I in 1:n
            for J in 1:I
                𝗚[I,J] += w*𝒒[I]*𝒒[J]
            end
        end
    end
    cholesky!(𝗚)
    inverse!(𝗚)
    UUᵀ!(𝗚)
    return 𝗚
end

function cal𝗚!(ap::ReproducingKernel{:Linear1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =  1.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Quadratic1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =  4.0/𝐿
    𝗚⁻¹[2] = -6.0/𝐿
    𝗚⁻¹[3] = 12.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Cubic1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =    9.0/𝐿
    𝗚⁻¹[2] =  -36.0/𝐿
    𝗚⁻¹[3] =  192.0/𝐿
    𝗚⁻¹[4] =   30.0/𝐿
    𝗚⁻¹[5] = -180.0/𝐿
    𝗚⁻¹[6] =  180.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Linear2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal∇𝗚₂!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝗚⁻¹∂ξ = get𝗚(ap,:∂∇̃²∂ξ)
    𝗚⁻¹∂η = get𝗚(ap,:∂∇̃²∂η)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴

    𝗚⁻¹∂ξ[1] =   9.0/𝐴
    𝗚⁻¹∂ξ[2] = -12.0/𝐴
    𝗚⁻¹∂ξ[3] =  24.0/𝐴
    𝗚⁻¹∂ξ[4] = -12.0/𝐴
    𝗚⁻¹∂ξ[5] =  12.0/𝐴
    𝗚⁻¹∂ξ[6] =  24.0/𝐴

    𝗚⁻¹∂η[1] =   9.0/𝐴
    𝗚⁻¹∂η[2] = -12.0/𝐴
    𝗚⁻¹∂η[3] =  24.0/𝐴
    𝗚⁻¹∂η[4] = -12.0/𝐴
    𝗚⁻¹∂η[5] =  12.0/𝐴
    𝗚⁻¹∂η[6] =  24.0/𝐴

    return 𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η
end

function cal𝗚₂!(ap::ReproducingKernel{:Quartic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end

function cal∇𝗚₂!(ap::ReproducingKernel{:Quartic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝗚⁻¹∂ξ = get𝗚(ap,:∂∇̃²∂ξ)
    𝗚⁻¹∂η = get𝗚(ap,:∂∇̃²∂η)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴

    𝗚⁻¹∂ξ[1] =   36.0/𝐴
    𝗚⁻¹∂ξ[2] = -120.0/𝐴
    𝗚⁻¹∂ξ[3] =  600.0/𝐴
    𝗚⁻¹∂ξ[4] = -120.0/𝐴
    𝗚⁻¹∂ξ[5] =  300.0/𝐴
    𝗚⁻¹∂ξ[6] =  600.0/𝐴
    𝗚⁻¹∂ξ[7] =   90.0/𝐴
    𝗚⁻¹∂ξ[8] = -540.0/𝐴
    𝗚⁻¹∂ξ[9] = -180.0/𝐴
    𝗚⁻¹∂ξ[10] =  540.0/𝐴
    𝗚⁻¹∂ξ[11] =  180.0/𝐴
    𝗚⁻¹∂ξ[12] = -720.0/𝐴
    𝗚⁻¹∂ξ[13] = -720.0/𝐴
    𝗚⁻¹∂ξ[14] =  540.0/𝐴
    𝗚⁻¹∂ξ[15] = 1440.0/𝐴
    𝗚⁻¹∂ξ[16] =   90.0/𝐴
    𝗚⁻¹∂ξ[17] = -180.0/𝐴
    𝗚⁻¹∂ξ[18] = -540.0/𝐴
    𝗚⁻¹∂ξ[19] =   90.0/𝐴
    𝗚⁻¹∂ξ[20] =  540.0/𝐴
    𝗚⁻¹∂ξ[21] =  540.0/𝐴

    𝗚⁻¹∂η[1] =   36.0/𝐴
    𝗚⁻¹∂η[2] = -120.0/𝐴
    𝗚⁻¹∂η[3] =  600.0/𝐴
    𝗚⁻¹∂η[4] = -120.0/𝐴
    𝗚⁻¹∂η[5] =  300.0/𝐴
    𝗚⁻¹∂η[6] =  600.0/𝐴
    𝗚⁻¹∂η[7] =   90.0/𝐴
    𝗚⁻¹∂η[8] = -540.0/𝐴
    𝗚⁻¹∂η[9] = -180.0/𝐴
    𝗚⁻¹∂η[10] =  540.0/𝐴
    𝗚⁻¹∂η[11] =  180.0/𝐴
    𝗚⁻¹∂η[12] = -720.0/𝐴
    𝗚⁻¹∂η[13] = -720.0/𝐴
    𝗚⁻¹∂η[14] =  540.0/𝐴
    𝗚⁻¹∂η[15] = 1440.0/𝐴
    𝗚⁻¹∂η[16] =   90.0/𝐴
    𝗚⁻¹∂η[17] = -180.0/𝐴
    𝗚⁻¹∂η[18] = -540.0/𝐴
    𝗚⁻¹∂η[19] =   90.0/𝐴
    𝗚⁻¹∂η[20] =  540.0/𝐴
    𝗚⁻¹∂η[21] =  540.0/𝐴

    return 𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η
end