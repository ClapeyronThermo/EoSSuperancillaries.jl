struct REFPROPSuperAnc
    vlsat::ChebyshevRangeV64
    vvsat::ChebyshevRangeV64
    psat::ChebyshevRangeV64
    Tc::Float64
    rhoc::Float64
    Rgas::Float64
    name::String
end

Base.show(io::IO,anc::REFPROPSuperAnc) = print(io,"REFPROPSuperAnc(\"$(anc.name)\")")
function Base.show(io::IO,::MIME"text/plain",anc::REFPROPSuperAnc)
    print(io,"REFPROPSuperAnc(\"$(anc.name)\")")
end

function normalisestring(str::String)
    res = Base.Unicode.normalize(str,casefold=true,stripmark=true)
    return replace(res, " " => "")
end

normalisestring(str) = normalisestring(String(str))

"""
    refprop_superanc(name::AbstractString)

Returns an instance of `REFPROPSuperAnc`.
It stores the ancillaries for the saturation pressure and densities for the input component.

"""
function refprop_superanc(name::AbstractString)
    file = REFPROP_SYNONYMS[normalisestring(name)]
    data = BSON.load(Base.@__DIR__()  * "/data/" * file)
    psat = ChebyshevRange(data[:jexpansions_p])
    vlsat = ChebyshevRange(data[:jexpansions_rhoL])
    vvsat = ChebyshevRange(data[:jexpansions_rhoV])
    
    meta = data[:meta]
    rhoc = meta[Symbol("rhocrittrue / mol/m^3")]
    Tc = meta[Symbol("Tcrittrue / K")]
    R = NaN
    return REFPROPSuperAnc(vlsat,vvsat,psat,Tc,rhoc,R,name)
end

"""
    p = ref_psat(fluid::REFPROPSuperanc, T)

Returns the saturation pressure for the reference fluid `fluid` at temperature `T`.

Inputs:
- `fluid`: Reference fluid (an instance of `REFPROPSuperAnc`)
- `T`: Temperature (Kelvin)

Outputs:
- `p` : Saturation Pressure `[Pa]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function ref_psat(anc::REFPROPSuperAnc, T)
    Tc = anc.Tc
    Θ = (Tc - T)/T
    return cheb_eval(anc.psat,T)
end

"""
    ρl = ref_rholsat(fluid::REFPROPSuperanc, T)

Returns the saturation liquid density for the reference fluid `fluid` at temperature `T`.

Inputs:
- `fluid`: Reference fluid (an instance of `REFPROPSuperAnc`)
- `T`: Temperature (Kelvin)

Outputs:
- `ρl` : Saturation Liquid Density `[mol/m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function ref_rholsat(anc::REFPROPSuperAnc, T)
    Tc = anc.Tc
    Θ = (Tc - T)/T
    return cheb_eval(anc.vlsat,T)
end

"""
    ρv = ref_rhovsat(fluid::REFPROPSuperanc, T)

Returns the saturation vapour density for the reference fluid `fluid` at temperature `T`.

Inputs:
- `fluid`: Reference fluid (an instance of `REFPROPSuperAnc`)
- `T`: Temperature (Kelvin)

Outputs:
- `ρv` : Saturation Vapour Density `[mol/m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function ref_rhovsat(anc::REFPROPSuperAnc, T)
    Tc = anc.Tc
    Θ = (Tc - T)/T
    return cheb_eval(anc.vvsat,T)
end


"""
    vl = ref_vlsat(fluid::REFPROPSuperanc, T)

Returns the saturation liquid volume for the reference fluid `fluid` at temperature `T`.

Inputs:
- `fluid`: Reference fluid (an instance of `REFPROPSuperAnc`)
- `T`: Temperature (Kelvin)

Outputs:
- `vl` : Saturation Liquid Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function ref_vlsat(anc::REFPROPSuperAnc, T)
    return 1/ref_rholsat(anc,T)
end

"""
    vv = ref_vvsat(fluid::REFPROPSuperanc, T)

Returns the saturation vapour volume for the reference fluid `fluid` at temperature `T`.

Inputs:
- `fluid`: Reference fluid (an instance of `REFPROPSuperAnc`)
- `T`: Temperature (Kelvin)

Outputs:
- `vv` : Saturation Vapour Volume `[m^3]`.  Returns `NaN` if the value is outside the range of the ancillary.
"""
function ref_vvsat(anc::REFPROPSuperAnc, T)
    return 1/ref_rhovsat(anc,T)
end
