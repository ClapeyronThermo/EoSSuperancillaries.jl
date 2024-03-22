# EoSSuperancillaries

Superancillary equations to calculate saturation pressures and critical points for selected equations of state.

This is a Julia port of [usnistgov/SAFTsuperanc](https://github.com/usnistgov/SAFTsuperanc)

## Introduction

Saturation calculation always involve iterative calculations. those calculations waste a lot of computing power. superancillaries provide direct equations that return some properties within `Float64` accuracy.

It can calculate the following properties:

- PC-SAFT: critical temperature, critical density, saturated densities.
- van der Wals: saturated densities, saturated pressures.
- Redlich-Kwong: saturated densities, saturated pressures.
- Peng-Robinson: saturated densities, saturated pressures.
- REFPROP reference equations: saturated densities, saturated pressures.

## Installation

You can install `EoSSuperancillaries` by running the following commands on the julia REPL:

```
julia> using Pkg

julia> Pkg.add("EoSSuperancillaries")
```

you can then load the package into your enviroment by writing in the same REPL:
```
julia> using EoSSuperancillaries
```
## Cubics

The most known cubic equations of state are defined by two parameters: `a` and `b`. The calculation of those parameters depends on each equation of state, values of the fluid at the critical point, and in the case of `a`, additional parameters that are used to define a temperature dependence. For example, if we take the (canonical) peng robinson equation of state:

```
a =  Ωa * (R * T)^2 / Pc * α(T)
b = Ωb * R * T / Pc
Ωa = 0.45723552892138218938
Ωb = 0.077796073903888455972
α = (1 + κ*(1 - √(T/Tc)))^2
κ = 0.37464 + 1.54226ω - 0.26992ω²
```

Where `Tc`,`Pc` is the critical temperature and pressure, and `ω` is the acentric factor. 

!!! note
    The choice of `α(T)` does not matter for cubic superancillaries, as we only use the `a` and `b` values. 

To use cubic superancillaries, the procedure is the following:
1. Calculate `a` and `b` values
2. pass those values to the adequate `EoSSuperancillaries` method.

```
using EoSSuperancillaries, Clapeyron
model = PR("water") #model for comparison.

#usage for cubics
a,b,_ = Clapeyron.cubic_ab(model,0.0,373.15)
#direct calculation using superancillaries, takes about 150 ns
p_eq = pr_psat(373.15,a,b) #96099.13581050883 
vl,vv = pr_vsat(373.15,a,b)
#iterative calculation, takes about 2 μs
p_eq2,vl2,vv2 = Clapeyron.saturation_pressure(model,373.15)[1] #96099.13581050835
```
All cubic superancillaries are implicitly defined in terms of `T̃ =RbT/a`, and have a range of application of `0.1T̃c < T̃ < T̃c`.

## PC-SAFT

The parturbed-chain statistical associating fluid theory (PC-SAFT) equation of state is defined, for non-associating compounds, in terms of the segment length (`m`), the monomer segment diameter (`σ`) and the reduced interaction energy (`ϵ`). Associating compounds require additional parameters (associating interaction energy and associating diameter, among with the number and type of sites), and are not supported by the current superancillary present in the package.

One PCSAFT superancillary calculation for the saturation volumes, involves the following:
1. calculate the critical temperature (we can do that with `pcsaft_tc(m,ϵ)`). in particular `T̃c = Tc/ϵ = pcsaft_tc(m,1.0)` is used. 
2. calculate the minimum temperature (`T̃min = T̃c*exp(-2.20078778)*m^0.37627892`)
3. Calculate `Θ = (T̃ - T̃min)/(T̃c - T̃min)`, where `T̃ = T/ϵ`
4. with `Θ` and `m`, calculate the reduced densities. `ρ̃l,ρ̃v`
5. reescale the densities: `ρ = ρ̃/(N_A*σ*σ*σ)` where `N_A` is the Avogadro's constant, equal to `6.02214076e23`

All those steps are done automatically by calling `pcsaft_vsat` (or `pcsaft_rhosat`), but some steps can be cached. For each compound, you can cache `T̃c` and `T̃min`, So you can do the following:

```julia
function my_pcsaft_rholsat(T,m,ϵ,σ)
    T̃ = T/ϵ
    Tc = pcsaft_tc(m,ϵ) #or from cached location
    T̃c = Tc/ϵ
    T̃min = T̃c*exp(-2.20078778)*m^0.37627892
    θ,status = pcsaft_theta(T̃,m,T̃c,T̃min)
    #you can also check isnan(θ) if you dont care about the exact out of range condition
    if status == :inrange
        ρ̃l = pcsaft_rholsat_reduced(θ,m)
        N_A = 6.02214076e23 #(or EoSSuperancillaries.N_A)
        return ρ̃l/(N_A*σ*σ*σ)
    end
```

This is useful in the cases where you want to benchmark your own saturation solver, or if you are working with extended precision arithmetic.

## REFPROP

Superancillaries for all 147 REFPROP 10.0 fluids have been fitted. those superancillaries cover from the critical point to the triple point of most fluids. for Helium, the minimum temperature is its lambda point. 

In this case, you need to instantiate a `REFPROPSuperAnc` object first and call the superancillaries with that object:

```julia
water = refprop_superanc("water")
p0 = ref_psat(water,373.15)
```

A table with synonyms accepted by `refprop_superanc` is found [here](@ref refprop_fluids)

## Notes on low reduced temperatures

With both cubics and PC-SAFT, the saturation curves are defined in a range `Tmin < T < Tc`. on low reduced temperatures, the magnitude of the proportion between volumes (`vv/vl`) reaches the limits of `Float64` arithmetic. for solving at temperatures below those ranges, it is recommended to use extended precision numbers (like `BigFloat`).

## References
1. Bell, I. H., & Deiters, U. K. (2021). Superancillary equations for cubic equations of state. Industrial & Engineering Chemistry Research, 60(27), 9983–9991. [doi:10.1021/acs.iecr.1c00847](https://doi.org/10.1021/acs.iecr.1c00847)
2. Bell, I. H., & Deiters, U. K. (2023). Superancillary equations for nonpolar pure fluids modeled with the PC-SAFT equation of state. Industrial & Engineering Chemistry Research. [doi:10.1021/acs.iecr.2c02916](https://doi.org/10.1021/acs.iecr.2c02916)
3. Bell, I. H. (2024). Superancillary equations for the multiparameter equations of state in REFPROP 10.0. Journal of Physical and Chemical Reference Data, 53(1). [doi:10.1063/5.0191228](https://doi.org/10.1063/5.0191228)