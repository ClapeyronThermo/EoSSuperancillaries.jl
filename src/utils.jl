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

function searchsortedfirst(xx::Tuple,x)
    for (i,xi) in pairs(xx)
        if xi >= x
            return i
        end
    end
    return length(xx)
end

searchsortedfirst(xx,x) = Base.searchsortedfirst(xx,x)

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

function cheb_xrange(x̃range,x̃,reverse = false)
    if !reverse
        imax = searchsortedfirst(x̃range, x̃)
        imin = imax - 1
        x̃minᵢ = x̃range[imin]
        x̃maxᵢ = x̃range[imax]
        i = imin
    else
        imax = searchsortedfirst(x̃range, x̃)
        imin = imax - 1
        x̃minᵢ = x̃range[imin]
        x̃maxᵢ = x̃range[imax]
        i = imin
    end
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