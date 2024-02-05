function _cubic_vsat(T,a,b,CV)
    T̃ = T*R̄*b/a
    ρ = chev_eval(CV,T̃)
    return b/ρ
end

function _cubic_vsat2(T,a,b,CV,VL)
    T̃ = T*R̄*b/a
    ρl = chev_eval(CL,T̃)
    ρv = chev_eval(CV,T̃)
    return b/ρl,b/ρv
end

function _cubic_psat(T,a,b,CP)
    T̃ = T*R̄*b/a
    p̃ = chev_eval(CP,T̃)
    return p̃*a/(b*b)    
end

"""
    vl = vdw_vlsat(T,a,b)

Returns the saturation liquid volume of the vdW cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function vdw_vlsat(T,a,b)
    return _cubic_vsat(T,a,b,vdW_vl)
end

"""
    vv = vdw_vvsat(T,a,b)

Returns the saturation vapour volume of the vdW cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vv` : Saturation Vapour Volume `[m^3]`. Returns `NaN` if the value is outside the range of the ancillary.
"""
function vdw_vvsat(T,a,b)
    return _cubic_vsat(T,a,b,vdW_vv)
end

"""
    p = vdw_psat(T,a,b)

Returns the saturation pressure of the vdW cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `p` : Saturation pressure `[Pa]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function vdw_psat(T,a,b)
    return _cubic_psat(T,a,b,vdW_p)
end


"""
    vl = rk_vlsat(T,a,b)

Returns the saturation liquid volume of the Redlich-Kwong cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function rk_vlsat(T,a,b)
    return _cubic_vsat(T,a,b,RK_vl)
end

"""
    vv = rk_vvsat(T,a,b)

Returns the saturation vapour volume of the Redlich-Kwong cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vv` : Saturation Vapour Volume `[m^3]`. Returns `NaN` if the value is outside the range of the ancillary.
"""
function rk_vvsat(T,a,b)
    return _cubic_vsat(T,a,b,RK_vv)
end

"""
    p = rk_psat(T,a,b)

Returns the saturation pressure of the Redlich-Kwong cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `p` : Saturation pressure `[Pa]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function rk_psat(T,a,b)
    return _cubic_psat(T,a,b,RK_p)
end

"""
    vl = pr_vlsat(T,a,b)

Returns the saturation liquid volume of the Peng-Robinson cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function pr_vlsat(T,a,b)
    return _cubic_vsat(T,a,b,PR_vl)
end

"""
    vv = pr_vvsat(T,a,b)

Returns the saturation vapour volume of the Peng-Robinson cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vv` : Saturation Vapour Volume `[m^3]`. Returns `NaN` if the value is outside the range of the ancillary.
"""
function pr_vvsat(T,a,b)
    return _cubic_vsat(T,a,b,PR_vv)
end

"""
    p = pr_psat(T,a,b)

Returns the saturation pressure of the Peng-Robinson cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `p` : Saturation pressure `[Pa]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function pr_psat(T,a,b)
    return _cubic_psat(T,a,b,PR_p)
end

