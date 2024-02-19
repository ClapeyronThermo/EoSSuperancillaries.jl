

function pcsaft_tc(m,ϵ)
    m⁻¹ = 1/m*oneunit(ϵ)*1.0
    T̃c = cheb_eval(PCSAFTsuperanc.Tc,m⁻¹)
    return T̃c*ϵ
end

function pcsaft_rhoc(m,ϵ,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(PCSAFTsuperanc.rho_c,m⁻¹)
    return ρ̃c/(N_A*σ*σ*σ)
end

function DCT_mat(x::Val{N}) where N
    v = @MMatrix zeros(N+1,N+1)
    for ii in 1:(N+1)
        for jj in 1:ii
            i,j = ii - 1, jj - 1
            p_i = (i == 0 || i == N) ? 2 : 1
            p_j = (j == 0 || j == N) ? 2 : 1
            vij = (2.0 / (N*p_i*p_j)) * cospi(jj*ii/N)
            v[ii,jj] = vij
            v[jj,ii] = vij
        end
    end
    return v
end

const V16 = DCT_mat(Val{16}())

function pcsaft_vsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    ρ̃l,ρ̃v = _pcsaft_rhosat(T̃,m)
    N_Aσ3 = N_A*σ*σ*σ
    vl = N_Aσ3/ρ̃l
    vv = N_Aσ3/ρ̃v
    return vl,vv
end

function _pcsaft_rhosat(T̃,m)
    _0 = zero(T̃+m+1.0)
    TYPE = typeof(_0)
    if !(1.0 <= m <= 64)
        nan = zero(TYPE)/zero(TYPE)
        return nan,nan
    end
    m⁻¹ = 1/m
    T̃c = cheb_eval(PCSAFTsuperanc.Tc,m⁻¹)
    T̃min = T̃c*exp(-2.20078778)*m^0.37627892
    if !(T̃min <= T̃ <= T̃c)
        nan = zero(TYPE)/zero(TYPE)
        return nan,nan
    end
    Θ = (T̃ - T̃min)/(T̃c - T̃min)
    vdata = PCSAFTsuperanc.WDomainComplete
    Wedges = vdata[1]
    WIntervals = vdata[2]
    w_idx = searchsortedfirst(Wedges,m⁻¹) - 1
    
    sat_data = WIntervals[w_idx]
    edges,liq,vap = sat_data

    ρ̃l_points = mvec(TYPE,Val{17}())
    ρ̃v_points = mvec(TYPE,Val{17}())

    for i in 1:17
        ρ̃l_points[i] = cheb_eval(liq[i],Θ)
        ρ̃v_points[i] = cheb_eval(vap[i],Θ)
    end   

    cheb_nodes_ρ̃l = V16*ρ̃l_points
    cheb_nodes_ρ̃v = V16*ρ̃v_points

    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,m⁻¹)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,m⁻¹)
    
    return ρ̃l,ρ̃v
end

function mvec(::Type{T},::Val{N}) where {T,N}
    if isbitstype(T) # MVector does not work on non bits types, like BigFloat
        return @MVector zeros(T,N)
    else
        return zeros(T,N)
    end
end
