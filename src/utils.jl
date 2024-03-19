#gas constant
const R̄ = 8.31446261815324

#Avogadro constant
const N_A = 6.02214076e23

#a struct that stores a vector of Chebyshev coefficients and a range where those coeff
#are applicable.
struct ChebyshevRange{R,T}
    range::R
    coeffs::T
end

const ChebishevRangeV64 = ChebyshevRange{Vector{Float64},Vector{Vector{Float64}}}

#overload of searchsortedfirst to also match the first element.
searchsortedfirst(xx::Tuple,x) = searchsortedfirst(SVector(xx),x)
searchsortedfirst(xx,x) = Base.searchsortedfirst(xx,x) + isequal(x,first(xx))

#evaluation of ranges of chebyshev coefficients
function cheb_eval(Base.@specialize(cheb), Base.@specialize(x̃))
    Cₙ = cheb.coeffs
    x̃range = cheb.range
    x̃min = first(x̃range)
    x̃max = last(x̃range)
    if !(x̃min <= x̃ <= x̃max)
        R = promote_type(eltype(x̃range), eltype(eltype(Cₙ)),typeof(x̃))
        #x is not in range
        return zero(R)/zero(R)
    end
    x̄,i = cheb_xrange(x̃range,x̃)
    Cₙi = Cₙ[i]
    return cheb_eval(Cₙi,x̄)
end

function cheb_xrange(x̃range,x̃)
        imax = searchsortedfirst(x̃range, x̃)
        imin = imax - 1
        x̃minᵢ = x̃range[imin]
        x̃maxᵢ = x̃range[imax]
        i = imin
    x̄ = (2*x̃ - (x̃maxᵢ + x̃minᵢ)) / (x̃maxᵢ - x̃minᵢ)
    return x̄,i
end

function cheb_eval(cheb::AbstractVector{T},x::S) where {T,S}
    R = promote_type(T, S)
    l = length(cheb)
    l == 0 && return zero(R)
    l == 1 && return R(cheb[1])
    c0 = cheb[l - 1]
    c1 = cheb[l]
    for i in (l-2):-1:1
        c0, c1 = cheb[i] - c1, c0 + c1 * 2x
    end
    return R(c0 + c1 * x)
end

function ChebyshevRange(data::Vector{Tuple{Vector{Float64},Float64,Float64}})
    n = length(data)
    c = Vector{Vector{Float64}}(undef,n)
    r = Vector{Float64}(undef,n + 1)
    r[1] = data[1][2] #first value of the range, the minimum
    for i in 1:n
        coeff,_,xmax = data[i]
        r[i+1] = xmax
        c[i] = coeff
    end
    return ChebyshevRange(r,c)
end

function ChebyshevRange(data::AbstractVector{T}) where T<:AbstractDict
    n = length(data)
    c = Vector{Vector{Float64}}(undef,n)
    r = Vector{Float64}(undef,n + 1)
    r[1] = data[1][:xmin] #first value of the range, the minimum
    for i in 1:n
        data_i = data[i]
        coeff = data_i[:coeff]
        xmax = data_i[:xmax]
        r[i+1] = xmax
        c[i] = coeff
    end
    return ChebyshevRange(r,c)
end

function dct_mat(::Val{N}) where N
    L = @MMatrix zeros(N+1,N+1)
    for j = 0:N
        for k = j:N
            p_j = (j == 0 || j == N) ? 2 : 1
            p_k = (k == 0 || k == N) ? 2 : 1
            L[j+1, k+1] = 2.0 / (p_j * p_k * N) * cos((j * pi * k) / N)
            L[k+1, j+1] = L[j+1, k+1]
        end
    end
    return L
end
