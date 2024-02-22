
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

pcsaft_theta(T̃,m) = _pcsaft_theta(T̃,1/m)
pcsaft_theta(T̃,m,T̃c) = _pcsaft_theta(T̃,1/m,T̃c)
pcsaft_theta(T̃,m,T̃c,T̃min) = _pcsaft_theta(T̃,1/m,T̃c,T̃min)

function _pcsaft_theta(T̃,m⁻¹,T̃c = cheb_eval(PCSAFTsuperanc.Tc,m⁻¹), T̃min = T̃c*exp(-2.20078778)*m⁻¹^-0.37627892)
    _0 = zero(T̃+m⁻¹+1.0)
    0.015625 <= m⁻¹ <= 1.0 || return _0/_0
    T̃min <= T̃ <= T̃c || return _0/_0
    Θ = (T̃ - T̃min)/(T̃c - T̃min)
    return Θ
end

function _pcsaft_rhosat(T̃,m)
    m⁻¹ = 1/m
    Θ = _pcsaft_theta(T̃,m⁻¹)
    isnan(Θ) && return Θ,Θ
    vdata = PCSAFTsuperanc.WDomainComplete
    Wedges = vdata[1]
    WIntervals = vdata[2]
    w_idx = searchsortedfirst(Wedges,m⁻¹) - 1
    sat_data = WIntervals[w_idx]
    edges,liq,vap = sat_data
    ρ̃l_points = svec17(liq,Θ)
    ρ̃v_points = svec17(vap,Θ)
    cheb_nodes_ρ̃l = Lmat*ρ̃l_points
    cheb_nodes_ρ̃v = Lmat*ρ̃v_points
    b,a = first(edges),last(edges)
    w = (2*m⁻¹ - (b + a))/(b - a)
    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,w)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,w)
    return ρ̃l,ρ̃v
end

function svec17(cheb,Θ)
    x1 = cheb_eval(cheb[1],Θ)
    x2 = cheb_eval(cheb[2],Θ)
    x3 = cheb_eval(cheb[3],Θ)
    x4 = cheb_eval(cheb[4],Θ)
    x5 = cheb_eval(cheb[5],Θ)
    x6 = cheb_eval(cheb[6],Θ)
    x7 = cheb_eval(cheb[7],Θ)
    x8 = cheb_eval(cheb[8],Θ)
    x9 = cheb_eval(cheb[9],Θ)
    x10 = cheb_eval(cheb[10],Θ)
    x11 = cheb_eval(cheb[11],Θ)
    x12 = cheb_eval(cheb[12],Θ)
    x13 = cheb_eval(cheb[13],Θ)
    x14 = cheb_eval(cheb[14],Θ)
    x15 = cheb_eval(cheb[15],Θ)
    x16 = cheb_eval(cheb[16],Θ)
    x17 = cheb_eval(cheb[17],Θ)
    return SVector((x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17))
end