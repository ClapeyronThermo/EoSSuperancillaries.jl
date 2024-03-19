function chebnorm(M,coef::AbstractVector)
    coef_hi = @view coef[1:M]
    coef_lo = @view coef[(end - M + 1):end]
    return norm(coef_lo)/norm(coef_hi)
end

function chebnorm(M,coef::Tuple)
   return chebnorm(M,SVector(coef))
end

chebnorm(M) = Base.Fix1(chebnorm,M)

function chebpoints(N)
    return [cos(i*pi/N) for i in 1:N]
end
"""
dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)

Given a function `f`, with limits `xmin` and `xmax`,returns a `ChebishevRange`, with a N-range Chebishev interpolation at each interval. 
The intervals are subdivided until `max_refine_passes` is reached or when the M-norm of chebishev coefficients is below `tol`.
"""
function dyadic_splitting(f,N,xmin,xmax;M = 3,tol = 1.0e-12,max_refine_passes = 16,threaded = false)
    return nothing
end