

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
            vij = 2.0 / (N*p_i*p_j) * cospi(jj*ii/N)
            v[ii,jj] = vij
            v[jj,ii] = vij
        end
    end
    return v
end

const mat17 = DCT_mat(Val{17}())
const V16 = DCT_mat(Val{16}())

function pcsaft_vsat(T,m,ϵ,σ)
    _0 = zero(T+m+ϵ+σ+1.0)
    TYPE = typeof(_0)
    T̃ = T/ϵ
    if !(1.0 <= m <= 64)
        nan = zero(TYPE)/zero(TYPE)
        return nan,nan
        #return StaticArrays.sacollect(SVector{17, TYPE}, nan for i in 1:17)
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
    #TODO: find the allocations
    static_edge = StaticArrays.sacollect(SVector{17, Float64},edges)
    ρ̃l_points = StaticArrays.sacollect(SVector{17, typeof(Θ)}, cheb_eval(liq[i],Θ) for i in 1:17)
    ρ̃v_points = StaticArrays.sacollect(SVector{17, typeof(Θ)}, cheb_eval(vap[i],Θ) for i in 1:17)
    
    cheb_nodes_ρ̃l = V16*ρ̃l_points
    cheb_nodes_ρ̃v = V16*ρ̃v_points
    #w,_ = cheb_xrange(reverse(edges),m⁻¹) #TODO: implement reverse
    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,m⁻¹)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,m⁻¹)

    ρl = ρ̃l/(N_A*σ*σ*σ)
    ρv = ρ̃v/(N_A*σ*σ*σ)
    return 1/ρl,1/ρv
    #static_edge = ntuple(i -> edges[i],Val{17}())
    #ρ̃l_points = NTuple{17, TYPE}(cheb_eval(liq[i],Θ) for i in 1:17) #ntuple(i -> cheb_eval(liq[i],Θ),Val{17}())
    #ρ̃v_points = NTuple{17, TYPE}(cheb_eval(vap[i],Θ) for i in 1:17)

    #return ρ̃l_m
    #ρ̃l_m = map(x -> cheb_eval(first(x),Θ), PCSAFT_vsat)
    #ρ̃v_m = map(x -> cheb_eval(last(x),Θ), PCSAFT_vsat)
end
