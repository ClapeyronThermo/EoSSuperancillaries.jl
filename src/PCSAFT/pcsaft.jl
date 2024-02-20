
"""
    Tc = pcsaft_tc(m,ϵ)

Returns the critical temperature of the PCSAFT equation of state.

Inputs:
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`

Outputs:
- `Tc` : Critical Temperature `[K]`.  Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).

## Examples
```julia-repl
julia> m, ϵ = 1.0, 150.03 #values for methane
(1.0, 150.03)

julia> Tc = pcsaft_tc(m, ϵ)
191.40058128833536
```
"""
function pcsaft_tc(m,ϵ)
    m⁻¹ = 1/m*oneunit(ϵ)*1.0
    T̃c = cheb_eval(PCSAFTsuperanc.Tc,m⁻¹)
    return T̃c*ϵ
end

"""
    rhoc = pcsaft_rhoc(m,ϵ,σ)

Returns the critical density of the PCSAFT equation of state.

Inputs:
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `rhoc` : Critical density `[mol/m^3]`. Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).
"""
function pcsaft_rhoc(m,ϵ,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(PCSAFTsuperanc.rho_c,m⁻¹)
    return ρ̃c/(N_A*σ*σ*σ)
end

"""
    vc = pcsaft_vc(m,ϵ,σ)

Returns the critical volume of the PCSAFT equation of state. Equal to `1/pcsaft_rhoc(m,ϵ,σ)`

Inputs:
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `vc` : Critical volume `[m^3/mol]`. Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).
"""
function pcsaft_vc(m,ϵ,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(PCSAFTsuperanc.rho_c,m⁻¹)
    return (N_A*σ*σ*σ)/ρ̃c
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

const Lmat = dct_mat(Val{16}())

"""
    vc = pcsaft_vsat(T,m,ϵ,σ)

Returns the saturation volume of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `vl` : saturation liquid volume `[m^3/mol]`. 
- `vv` : saturation vapour volume `[m^3/mol]`.

Returns `NaN,NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).
"""
function pcsaft_vsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    ρ̃l,ρ̃v = _pcsaft_rhosat(T̃,m)
    N_Aσ3 = N_A*σ*σ*σ
    vl = N_Aσ3/ρ̃l
    vv = N_Aσ3/ρ̃v
    return vl,vv
end

"""
    vc = pcsaft_rhosat(T,m,ϵ,σ)

Returns the saturation volumes of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `rhol` : saturation liquid volume `[m^3/mol]`. 
- `rhov` : saturation vapour volume `[m^3/mol]`.

Returns `NaN,NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).
"""
function pcsaft_rhosat(T,m,ϵ,σ)
    T̃ = T/ϵ
    ρ̃l,ρ̃v = _pcsaft_rhosat(T̃,m)
    N_Aσ3 = N_A*σ*σ*σ
    ρl = ρ̃l/N_Aσ3
    ρv = ρ̃v/N_Aσ3
    return ρl,ρv
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
    w_idx += isequal(m,64) #fix from when the value is exactly 64
    sat_data = WIntervals[w_idx]
    edges,liq,vap = sat_data

    ρ̃l_points = mvec(TYPE,Val{17}())
    ρ̃v_points = mvec(TYPE,Val{17}())

    for i in 1:17
        ρ̃l_points[i] = cheb_eval(liq[i],Θ)
        ρ̃v_points[i] = cheb_eval(vap[i],Θ)
    end   

    cheb_nodes_ρ̃l = Lmat*ρ̃l_points
    cheb_nodes_ρ̃v = Lmat*ρ̃v_points
    b,a = first(edges),last(edges)
    w = (2*m⁻¹ - (b + a))/(b - a)
    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,w)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,w)
    return ρ̃l,ρ̃v
end

function mvec(::Type{T},::Val{N}) where {T,N}
    if isbitstype(T) # MVector does not work on non bits types, like BigFloat
        return @MVector zeros(T,N)
    else
        return zeros(T,N)
    end
end
