
[![Build Status](https://github.com/ClapeyronThermo/EoSSuperancillaries.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ClapeyronThermo/EoSSuperancillaries.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://clapeyronthermo.github.io/EoSSuperancillaries.jl/dev)
[![project chat](https://img.shields.io/badge/zulip-join_chat-brightgreen.svg)](https://julialang.zulipchat.com/#narrow/stream/265161-Clapeyron.2Ejl)

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

## Usage

You can install `EoSSuperancillaries` by running the following commands on the julia REPL:

```
using Pkg
Pkg.add(url="https://github.com/ClapeyronThermo/EoSSuperancillaries.jl")
```

```
using EoSSuperancillaries, Clapeyron
model = PR("water") #model for comparison.

#usage for cubics
a,b,_ = Clapeyron.cubic_ab(model,0.0,373.15)
#direct calculation using superancillaries, takes about 150 ns
p_eq = pr_psat(373.15,a,b) #96099.13581050883 
vl,vv = pr_vlsat(373.15,a,b), pr_vvsat(373.15,a,b)
#iterative calculation, takes about 2 μs
p_eq2,vl2,vv2 = Clapeyron.saturation_pressure(model,373.15)[1] #96099.13581050835

#usage for PC-SAFT
saft = PCSAFT("propane")
m = saft.params.segment[1]
ϵ = saft.params.epsilon[1]
Tc = pcsaft_tc(m,ϵ) #150 ns
Tc2 = Clapeyron.crit_pure(saft) #around 100 μs
```

## Notes
This is not a library to calculate general properties of equations of state. in particular, you would need a PC-SAFT pressure implementation for calculation of saturation pressures and critical pressure.

While cubics can accept any alpha function for the calculation of `a`, PCSAFT superancillaries are more restricted, in the sense that they are fitted for components without association terms.

Superancillaries have problems on really low reduced temperatures, but those problems are caused by the artifacts of `Float64` arithmetic. PCSAFT in particular has non-physical behaviour at low reduced temperatures, and at high reduced temperatures in conjuction with a high value of the segment length.

## References
1. Bell, I. H., & Deiters, U. K. (2021). Superancillary equations for cubic equations of state. Industrial & Engineering Chemistry Research, 60(27), 9983–9991. [doi:10.1021/acs.iecr.1c00847](https://doi.org/10.1021/acs.iecr.1c00847)
2. Bell, I. H., & Deiters, U. K. (2023). Superancillary equations for nonpolar pure fluids modeled with the PC-SAFT equation of state. Industrial & Engineering Chemistry Research. [doi:10.1021/acs.iecr.2c02916](https://doi.org/10.1021/acs.iecr.2c02916)