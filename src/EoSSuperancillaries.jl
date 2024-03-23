module EoSSuperancillaries
using StaticArrays
using LinearAlgebra
using BSON
include("utils.jl")
include("dyadic_splitting.jl")
include("Cubic/consts.jl")
include("Cubic/cubic.jl")
include("PCSAFT/consts.jl")
include("PCSAFT/pcsaft.jl")
include("REFPROP/synonyms.jl")
include("REFPROP/refprop.jl")

export pr_psat,pr_vlsat,pr_vvsat,pr_vsat
export vdw_psat,vdw_vlsat,vdw_vvsat,vdw_vsat
export rk_psat,rk_vlsat,rk_vvsat,rk_vsat

export pcsaft_tc,pcsaft_rhoc,pcsaft_vc,pcsaft_acentric,pcsaft_m_from_acentric
export pcsaft_theta
export pcsaft_rholsat_reduced,pcsaft_rhovsat_reduced,pcsaft_rhosat_reduced
export pcsaft_vlsat,pcsaft_vvsat,pcsaft_vsat
export pcsaft_rhosat

export refprop_superanc
export ref_vlsat,ref_vvsat,ref_psat,ref_rholsat,ref_rhovsat

end #module
