function _cubic_vsat(T,a,b,CV)
    T̃ = T*R̄*b/a
    ρ = cheb_eval(CV,T̃)
    return b/ρ
end

function _cubic_vsat2(T,a,b,CV,VL)
    T̃ = T*R̄*b/a
    ρl = cheb_eval(CL,T̃)
    ρv = cheb_eval(CV,T̃)
    return b/ρl,b/ρv
end

function _cubic_psat(T,a,b,CP)
    T̃ = T*R̄*b/a
    p̃ = cheb_eval(CP,T̃)
    return p̃*a/(b*b)    
end

_vdw_p() = vdW_p
_vdw_vl() = vdW_vl
_vdw_vv() = vdW_vv

_rk_p() = RK_p
_rk_vl() = RK_vl
_rk_vv() = RK_vv

_pr_p() = PR_p
_pr_vl() = PR_vl
_pr_vv() = PR_vv
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
    return _cubic_vsat(T,a,b,_vdw_vl())
end

"""
    vl,vv = vdw_vsat(T,a,b)

Returns the saturation liquid and vapour volumes of the vdW cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
- `vv` : Saturation Vapour Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function vdw_vsat(T,a,b)
    return _cubic_vsat2(T,a,b,_vdw_vl(),_vdw_vv())
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
    return _cubic_vsat(T,a,b,_vdw_vv())
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
    return _cubic_psat(T,a,b,_vdw_p())
end

"""
    vl,vv = rk_vsat(T,a,b)

Returns the saturation liquid and vapour volumes of the Redlich-Kwong cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
- `vv` : Saturation Vapour Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function rk_vsat(T,a,b)
    return _cubic_vsat2(T,a,b,_rk_vl(),_rk_vv())
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
    return _cubic_vsat(T,a,b,_rk_vl())
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
    return _cubic_vsat(T,a,b,_rk_vv())
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
    return _cubic_psat(T,a,b,_rk_p())
end

"""
    vl,vv = pr_vsat(T,a,b)

Returns the saturation liquid and vapour volumes of the Peng-Robinson cubic equation of state.

Inputs:
- `T`: Temperature (Kelvin)
- `a`: Atraction parameter `[m^6*Pa/mol]`
- `b`: Covolume `[m^3/mol]`

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
- `vv` : Saturation Vapour Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function pr_vsat(T,a,b)
    return _cubic_vsat2(T,a,b,_pr_vl(),_pr_vv())
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
    return _cubic_vsat(T,a,b,_pr_vl())
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
    return _cubic_vsat(T,a,b,_pr_vv())
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
    return _cubic_psat(T,a,b,_pr_p())
end

