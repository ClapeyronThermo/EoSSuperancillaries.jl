function chebnorm(M,coef::AbstractVector)
    coef_hi = @view coef[1:M]
    coef_lo = @view coef[(end - M + 1):end]
    return norm(coef_lo)/norm(coef_hi)
end

const CHEB_MAT_STORE = Dict{Tuple{Int,DataType},Matrix}()
const CHEB_POINT_STORE = Dict{Tuple{Int,DataType},Vector}()

function chebnorm(M,coef::Tuple)
   return chebnorm(M,SVector(coef))
end

chebpoints(N) = chebpoints(N,Float64)

function chebpoints(N,::Type{T}) where T
    n = T(N)
    return [cos(i*(pi/n)) for i in 0:N]
end

chebnorm(M) = Base.Fix1(chebnorm,M)

struct ChebTransformed{F,T}
    f::F
    xmin::T
    xmax::T
end

function (f̄::ChebTransformed{F,T})(x̄) where {F,T}
    x̄max,x̄min = f̄.xmax,f̄.xmin
    x = 0.5*x̄*(x̄max - x̄min) + 0.5*(x̄max + x̄min)
    convert(T,f̄.f(x))
end

function cheb_interp(f,xmin::T,xmax::T,N::Int) where T
    f̄ = ChebTransformed(f,xmin,xmax)
    key = (N,typeof(xmin))
    x̄ = get!(CHEB_POINT_STORE,key) do
        chebpoints(N,T)
    end
    M = get(CHEB_MAT_STORE,key) do
        dct_mat(N,T)
    end
    coef = M*map(f̄,x̄)
    return ChebyshevRange([xmin,xmax],[coef])
end

"""
dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)

Given a function `f`, with limits `xmin` and `xmax`,returns a `ChebishevRange`, with a N-range Chebishev interpolation at each interval.
The intervals are subdivided until `max_refine_passes` is reached or when the M-norm of chebishev coefficients is below `tol`.
"""
function dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)
    #try to use 1 Chebyshev interpolation
    cheb = cheb_interp(f,xmin,xmax,N)
    if chebnorm(M,cheb.coeffs[1]) < tol
        return cheb
    end
    if max_refine_passes == 0
        @warn "max refine passes reached at x ∈ [$xmin,$xmax]"
        return cheb
    end
    #split domain
    xmid = 0.5*(xmin + xmax)
    left = dyadic_splitting(f,N,xmin,xmid,M = M,tol = tol,max_refine_passes = max_refine_passes - 1,threaded = threaded)
    right = dyadic_splitting(f,N,xmid,xmax,M = M,tol = tol,max_refine_passes = max_refine_passes - 1,threaded = threaded)
    return merge_chebs!(left,right)
end