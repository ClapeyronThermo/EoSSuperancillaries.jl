const Lmat16 = dct_mat(Val{16}())

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
    ω = pcsaft_acentric(m)

Returns the acentric factor of the PCSAFT equation of state.

Inputs:
- `m`: Segment length (no units)

Outputs:
- `ω` : acentric factor (no units).  Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).

## Examples
```julia-repl
julia> m = 1.0
(1.0, 150.03)

julia> Tc = pcsaft_acentric(m)
191.40058128833536
```
"""
function pcsaft_acentric(m)
    m⁻¹ = 1/m
    return cheb_eval(PCSAFTsuperanc.acentric,m⁻¹)
end

"""
    m = pcsaft_m_from_acentric(w)

Given the acentric factor of a component, returns the segment length that corresponds to the real acentric factor calculated by the PCSAFT equation of state.

Inputs:
- `ω` : acentric factor (no units)


Outputs:
- `m`: Segment length (no units).  Returns `NaN` if the value is outside the range of the ancillary (0.005073686089814064 ≤ w ≤ 6.0465684112508296).
"""
function pcsaft_m_from_acentric(ω)
    return cheb_eval(PCSAFTsuperanc.m_from_acentric, ω)
end

"""
    rhoc = pcsaft_rhoc(m,σ)

Returns the critical density of the PCSAFT equation of state.

Inputs:
- `m`: Segment length (no units)
- `σ`: Monomer diameter `[m]`

Outputs:
- `rhoc` : Critical density `[mol/m^3]`. Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).
"""
function pcsaft_rhoc(m,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(PCSAFTsuperanc.rhoc,m⁻¹)
    return ρ̃c/(N_A*σ*σ*σ)
end

"""
    vc = pcsaft_vc(m,σ)

Returns the critical volume of the PCSAFT equation of state. Equal to `1/pcsaft_rhoc(m,ϵ,σ)`

Inputs:
- `m`: Segment length (no units)
- `σ`: Monomer diameter `[m]`

Outputs:
- `vc` : Critical volume `[m^3/mol]`. Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64).
"""
function pcsaft_vc(m,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(PCSAFTsuperanc.rhoc,m⁻¹)
    return (N_A*σ*σ*σ)/ρ̃c
end


"""
    vl,vv = pcsaft_vsat(T,m,ϵ,σ)

Returns the saturation volumes of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `vl` : saturation liquid volume `[m^3/mol]`.
- `vv` : saturation vapour volume `[m^3/mol]`.

Returns `NaN,NaN` if the value is outside the range of the ancillary (1 ≤ `m` ≤ 64 and `Tmin` < `T` < `Tc`).
"""
function pcsaft_vsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    θ,_ = _pcsaft_theta(T̃,m⁻¹)
    ρ̃l,ρ̃v = _pcsaft_rhosat(θ,m⁻¹)
    N_Aσ3 = N_A*σ*σ*σ
    vl = N_Aσ3/ρ̃l
    vv = N_Aσ3/ρ̃v
    return vl,vv
end

"""
    rhol,rhov = pcsaft_rhosat(T,m,ϵ,σ)

Returns the saturation densities of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `rhol` : saturation liquid density `[mol/m^3]`.
- `rhov` : saturation vapour density `[mol/m^3]`.

Returns `NaN,NaN` if the value is outside the range of the ancillary (1 ≤ `m` ≤ 64 and `Tmin` < `T` < `Tc`).
"""
function pcsaft_rhosat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    θ,_ = _pcsaft_theta(T̃,m⁻¹)
    ρ̃l,ρ̃v = _pcsaft_rhosat(θ,m⁻¹)
    N_Aσ3 = N_A*σ*σ*σ
    ρl = ρ̃l/N_Aσ3
    ρv = ρ̃v/N_Aσ3
    return ρl,ρv
end

"""
    vl = pcsaft_vlsat(T,m,ϵ,σ)

Returns the saturation liquid volume of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `vl` : saturation liquid volume `[m^3/mol]`.

Returns `NaN` if the value is outside the range of the ancillary (1 ≤ `m` ≤ 64 and `Tmin` < `T` < `Tc`).
"""
function pcsaft_vvsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    θ,_ = _pcsaft_theta(T̃,m⁻¹)
    ρ̃v = _pcsaft_rhovsat(θ,m⁻¹)
    N_Aσ3 = N_A*σ*σ*σ
    vv = N_Aσ3/ρ̃v
    return vv
end

"""
    vv = pcsaft_vvsat(T,m,ϵ,σ)

Returns the saturation vapour volume of the PCSAFT equation of state at the input temperature `T`.

Inputs:
- `T`: Saturation temperature (Kelvin)
- `m`: Segment length (no units)
- `ϵ`: Reduced interaction energy `[K]`
- `σ`: Monomer diameter `[m]`

Outputs:
- `vl` : saturation liquid volume `[m^3/mol]`.

Returns `NaN` if the value is outside the range of the ancillary (1 ≤ `m` ≤ 64 and `Tmin` < `T` < `Tc`).
"""
function pcsaft_vlsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    θ,_ = _pcsaft_theta(T̃,m⁻¹)
    ρ̃l = _pcsaft_rholsat(θ,m⁻¹)
    N_Aσ3 = N_A*σ*σ*σ
    vl = N_Aσ3/ρ̃l
    return vl
end

"""
    Θ,status = pcsaft_theta(T̃,m)
    Θ,status = pcsaft_theta(T̃,m,T̃c)

Calculates reduced scaled temperature parameter used in PCSAFT saturation volume superancillary

Inputs:
- `T̃`: Temperature divided by reduced interaction energy (`T/ϵ`) (no units)
- `m`: Segment length (no units)
- `T̃c`: (Optional) Reduced critical temperature (`Tc/ϵ`) (no units)

Outputs:
- `Θ` : scaled temperature parameter (no units)
- `status` : `Symbol` used to signal if the combination of `T̃`,`m` is valid or not. It can return one of the following:

    - `:inrange` : if the combination is valid
    - `:nonfinite` : if any input is not finite
    - `:below_mmin` : if `m` < 1.0
    - `:over_mmax` : if `m` > 64.0
    - `:below_Tmin` : if `T̃` < `T̃min`, where `T̃min` = `T̃c*exp(-2.20078778)*m^0.37627892`
    - `:over_Tmax` : if `T̃` > `T̃c`
"""
function pcsaft_theta end

pcsaft_theta(T̃,m) = _pcsaft_theta(T̃,1/m)
pcsaft_theta(T̃,m,T̃c) = _pcsaft_theta(T̃,1/m,T̃c)
pcsaft_theta(T̃,m,T̃c,T̃min) = _pcsaft_theta(T̃,1/m,T̃c,T̃min)

function _pcsaft_theta(T̃,m⁻¹,T̃c = cheb_eval(PCSAFTsuperanc.Tc,m⁻¹), T̃min = T̃c*exp(-2.20078778)*m⁻¹^-0.37627892)
    _0 = zero(T̃+m⁻¹+1.0)
    if isnan(T̃) || isnan(m⁻¹) || isnan(T̃c) || isnan(T̃min)
        return _0/_0,:nonfinite
    elseif 0.015625 > m⁻¹
        return _0/_0,:over_mmax
    elseif m⁻¹ > 1.0
        return _0/_0,:under_mmin
    elseif T̃min > T̃
        return _0/_0,:below_Tmin
    elseif T̃ > T̃c
        return _0/_0,:over_Tmax
    else
        Θ = (T̃ - T̃min)/(T̃c - T̃min)
        return Θ,:inrange
    end
end

"""
    ρ̃l = pcsaft_rholsat_reduced(Θ,m)

Calculates the reduced liquid saturation density. The actual density is `ρ = ρ̃/(N_A*σ^3)`

Inputs:

- `Θ`: reduced temperature parameter (no units)
- `m`: Segment length (no units)

Outputs:
- `ρ̃l` : reduced saturation liquid density (no units)

Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).
"""
function pcsaft_rholsat_reduced(Θ,m)
    m⁻¹ = 1/m
    return _pcsaft_rholsat(Θ,m⁻¹)
end

"""
    ρ̃v = pcsaft_rhovsat_reduced(Θ,m)

Calculates the reduced vapour saturation density. The actual density is `ρ = ρ̃/(N_A*σ^3)`

Inputs:

- `Θ`: reduced temperature parameter (no units)
- `m`: Segment length (no units)

Outputs:
- `ρ̃v` : reduced saturation vapour density (no units)

Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).
"""
function pcsaft_rhovsat_reduced(Θ,m)
    m⁻¹ = 1/m
    return _pcsaft_rhovsat(Θ,m⁻¹)
end

"""
    ρ̃l,ρ̃v = pcsaft_rhosat_reduced(Θ,m)

Calculates the reduced saturation densities. The actual densities are calculated as `ρᵢ = ρ̃ᵢ/(N_A*σ^3)`

Inputs:

- `Θ`: reduced temperature parameter (no units)
- `m`: Segment length (no units)

Outputs:
- `ρ̃l` : reduced saturation liquid density (no units)
- `ρ̃v` : reduced saturation vapour density (no units)

Returns `NaN` if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).
"""
function pcsaft_rhosat_reduced(Θ,m)
    m⁻¹ = 1/m
    return _pcsaft_rhosat(Θ,m⁻¹)
end

function _pcsaft_rhosat(Θ,m⁻¹)
    isnan(Θ) && return Θ,Θ
    vdata = PCSAFTsuperanc.WDomainComplete
    Wedges = vdata[1]
    WIntervals = vdata[2]
    w_idx = searchsortedfirst(Wedges,m⁻¹) - 1
    sat_data = WIntervals[w_idx]
    edges,liq,vap = sat_data
    ρ̃l_points = svec17(liq,Θ)
    ρ̃v_points = svec17(vap,Θ)
    cheb_nodes_ρ̃l = Lmat16*ρ̃l_points
    cheb_nodes_ρ̃v = Lmat16*ρ̃v_points
    b,a = first(edges),last(edges)
    w = (2*m⁻¹ - (b + a))/(b - a)
    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,w)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,w)
    return ρ̃l,ρ̃v
end

function _pcsaft_rholsat(Θ,m⁻¹)
    isnan(Θ) && return Θ
    vdata = PCSAFTsuperanc.WDomainComplete
    Wedges = vdata[1]
    WIntervals = vdata[2]
    w_idx = searchsortedfirst(Wedges,m⁻¹) - 1
    sat_data = WIntervals[w_idx]
    edges,liq,_ = sat_data
    ρ̃l_points = svec17(liq,Θ)
    cheb_nodes_ρ̃l = Lmat16*ρ̃l_points
    b,a = first(edges),last(edges)
    w = (2*m⁻¹ - (b + a))/(b - a)
    ρ̃l = cheb_eval(cheb_nodes_ρ̃l,w)
    return ρ̃l
end

function _pcsaft_rhovsat(Θ,m⁻¹)
    isnan(Θ) && return Θ
    vdata = PCSAFTsuperanc.WDomainComplete
    Wedges = vdata[1]
    WIntervals = vdata[2]
    w_idx = searchsortedfirst(Wedges,m⁻¹) - 1
    sat_data = WIntervals[w_idx]
    edges,_,vap = sat_data
    ρ̃v_points = svec17(vap,Θ)
    cheb_nodes_ρ̃v = Lmat16*ρ̃l_points
    b,a = first(edges),last(edges)
    w = (2*m⁻¹ - (b + a))/(b - a)
    ρ̃v = cheb_eval(cheb_nodes_ρ̃v,w)
    return ρ̃v
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
