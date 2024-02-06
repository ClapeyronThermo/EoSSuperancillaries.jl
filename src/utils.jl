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

struct ChebyshevBiRange{R,T}
    range::R
    coeffs1::T
    coeffs2::T
end

# a lot of chevishev interpolations share the same range.
Base.first(r::ChebyshevBiRange) = ChebyshevRange(r.range,r.coeffs1)
Base.last(r::ChebyshevBiRange) = ChebyshevRange(r.range,r.coeffs2)

#evaluation of ranges of chebyshev coefficients
function cheb_eval(Base.@specialize(cheb), Base.@specialize(x̃))
    Cₙ = cheb.coeffs
    x̃range = cheb.range
    x̃min = first(x̃range)
    x̃max = last(x̃range)
    if !(x̃min <= x̃ <= x̃max)
        #x is not in range
        nan = 0.0*zero(x̃)/zero(x̃)
        return nan
    end
    imax = searchsortedfirst(x̃range, x̃)
    imin = imax - 1
    x̃minᵢ = x̃range[imin]
    x̃maxᵢ = x̃range[imax]
    Cₙi = Cₙ[imin]
    x̄ = (2*x̃ - (x̃maxᵢ + x̃minᵢ)) / (x̃maxᵢ - x̃minᵢ)
    return cheb_eval(Cₙi,x̄)
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