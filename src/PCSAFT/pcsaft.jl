function pcsaft_tc(m,Base.@specialize(ϵ))
    m⁻¹ = 1/m*oneunit(ϵ)*1.0
    T̃c = cheb_eval(first(PCSAFT_critical),m⁻¹)::typeof(m⁻¹)
    return T̃c*ϵ
end

function pcsaft_rhoc(m,ϵ,σ)
    m⁻¹ = 1/m
    ρ̃c = cheb_eval(last(PCSAFT_critical),m⁻¹)
    return ρ̃c/(N_A*σ*σ*σ)
end

function pcsaft_vsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    T̃c = cheb_eval(first(PCSAFT_critical),m⁻¹)
    T̃min = T̃c*exp(-2.20078778)*m^0.37627892
    if !(T̃min <= T̃ <= T̃c)
        nan = zero(T̃)/zero(T̃)
        return nan,nan
    end
    Θ = (T̃ - T̃min)/(T̃c - T̃min)
    #ρ̃l_m = map(x -> cheb_eval(first(x),Θ), PCSAFT_vsat)
    #ρ̃v_m = map(x -> cheb_eval(last(x),Θ), PCSAFT_vsat)
    #TODO: implement interpolation here.
    Nσ3 = N_A*σ*σ*σ
    ρ̃l = cheb_eval(PCSAFT_vsat_m01 |> first,Θ)
    vl = Nσ3/ρ̃l
    ρ̃v = cheb_eval(PCSAFT_vsat_m01 |> last,Θ)
    vv = Nσ3/ρ̃v
    return vl,vv
end
