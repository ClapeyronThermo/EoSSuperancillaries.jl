function pcsaft_tc(m,Base.@specialize(ϵ))
    m⁻¹ = 1/m*oneunit(ϵ)*1.0
    T̃c = chev_eval(first(PCSAFT_critical),m⁻¹)::typeof(m⁻¹)
    return T̃c*ϵ
end

function pcsaft_rhoc(m,ϵ,σ)
    m⁻¹ = 1/m
    ρ̃c = chev_eval(last(PCSAFT_critical),m⁻¹)
    return ρ̃c/(N_A*σ*σ*σ)
end

function pcsaft_vsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    m⁻¹ = 1/m
    T̃c = chev_eval(first(PCSAFT_critical),m⁻¹)
    T̃min = T̃c*exp(-2.20078778)*m^0.37627892
    if !(T̃min <= T̃ <= T̃c)
        nan = zero(T̃)/zero(T̃)
        return nan,nan
    end
    Θ = (T̃ - T̃min)/(T̃c - T̃min)
    ρ̃l_m = map(x -> chev_eval(first(x),Θ), PCSAFT_vsat)
    ρ̃v_m = map(x -> chev_eval(last(x),Θ), PCSAFT_vsat)
    #TODO: implement interpolation here.
    Nσ3 = N_A*σ*σ*σ
    ρ̃l = chev_eval(PCSAFT_vsat_m01 |> first,Θ)
    vl = Nσ3/ρ̃l
    return vl
end

#=
f(x) = ((x-0.5)^3)/sqrt(x + 8)
y = x + 1
f2(z) = -((1/z-1.5)^3)/sqrt(1/z + 7)*exp(z)/z/z #0 to 2

f3(z) = -exp(3*log(1/z - 1.5) + 1)/sqrt(1/z + 7)/z/z

z = 1/y
dz = -1/y^2dy

inf to 0.5



x = 0:0.01:1
y = f.(x

=#

#[0.01947021484375, 0.019433272481244004, 0.019323865067999605, 0.019146197069551403, 0.0189070961674424, 0.01861575087525131, 0.01828335742916872, 0.017922689522919193, 0.017547607421875, 0.017172525320830807, 0.01681185741458128, 0.01647946396849869, 0.0161881186763076, 0.015949017774198597, 0.015771349775750395, 0.015661942362505996, 0.015625]