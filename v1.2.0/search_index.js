var documenterSearchIndex = {"docs":
[{"location":"rk/","page":"Redlich-Kwong","title":"Redlich-Kwong","text":"CurrentModule = EoSSuperancillaries","category":"page"},{"location":"rk/#Contents","page":"Redlich-Kwong","title":"Contents","text":"","category":"section"},{"location":"rk/","page":"Redlich-Kwong","title":"Redlich-Kwong","text":"Pages = [\"rk.md\"]","category":"page"},{"location":"rk/#Index","page":"Redlich-Kwong","title":"Index","text":"","category":"section"},{"location":"rk/","page":"Redlich-Kwong","title":"Redlich-Kwong","text":"Pages = [\"rk.md\"]","category":"page"},{"location":"rk/#Redlich-Kwong-superancillaries","page":"Redlich-Kwong","title":"Redlich-Kwong superancillaries","text":"","category":"section"},{"location":"rk/","page":"Redlich-Kwong","title":"Redlich-Kwong","text":"EoSSuperancillaries.rk_vsat\nEoSSuperancillaries.rk_vlsat\nEoSSuperancillaries.rk_vvsat\nEoSSuperancillaries.rk_psat","category":"page"},{"location":"rk/#EoSSuperancillaries.rk_vsat","page":"Redlich-Kwong","title":"EoSSuperancillaries.rk_vsat","text":"vl,vv = rk_vsat(T,a,b)\n\nReturns the saturation liquid and vapour volumes of the Redlich-Kwong cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\nvv : Saturation Vapour Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"rk/#EoSSuperancillaries.rk_vlsat","page":"Redlich-Kwong","title":"EoSSuperancillaries.rk_vlsat","text":"vl = rk_vlsat(T,a,b)\n\nReturns the saturation liquid volume of the Redlich-Kwong cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"rk/#EoSSuperancillaries.rk_vvsat","page":"Redlich-Kwong","title":"EoSSuperancillaries.rk_vvsat","text":"vv = rk_vvsat(T,a,b)\n\nReturns the saturation vapour volume of the Redlich-Kwong cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvv : Saturation Vapour Volume [m^3]. Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"rk/#EoSSuperancillaries.rk_psat","page":"Redlich-Kwong","title":"EoSSuperancillaries.rk_psat","text":"p = rk_psat(T,a,b)\n\nReturns the saturation pressure of the Redlich-Kwong cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\np : Saturation pressure [Pa].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/","page":"PC-SAFT","title":"PC-SAFT","text":"CurrentModule = EoSSuperancillaries","category":"page"},{"location":"pcsaft/#Contents","page":"PC-SAFT","title":"Contents","text":"","category":"section"},{"location":"pcsaft/","page":"PC-SAFT","title":"PC-SAFT","text":"Pages = [\"pcsaft.md\"]","category":"page"},{"location":"pcsaft/#Index","page":"PC-SAFT","title":"Index","text":"","category":"section"},{"location":"pcsaft/","page":"PC-SAFT","title":"PC-SAFT","text":"Pages = [\"pcsaft.md\"]","category":"page"},{"location":"pcsaft/#PC-SAFT-critical-point","page":"PC-SAFT","title":"PC-SAFT critical point","text":"","category":"section"},{"location":"pcsaft/","page":"PC-SAFT","title":"PC-SAFT","text":"EoSSuperancillaries.pcsaft_tc\nEoSSuperancillaries.pcsaft_rhoc\nEoSSuperancillaries.pcsaft_vc\n","category":"page"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_tc","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_tc","text":"Tc = pcsaft_tc(m,ϵ)\n\nReturns the critical temperature of the PCSAFT equation of state.\n\nInputs:\n\nm: Segment length (no units)\nϵ: Reduced interaction energy [K]\n\nOutputs:\n\nTc : Critical Temperature [K].  Returns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64).\n\nExamples\n\njulia> m, ϵ = 1.0, 150.03 #values for methane\n(1.0, 150.03)\n\njulia> Tc = pcsaft_tc(m, ϵ)\n191.40058128833536\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_rhoc","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_rhoc","text":"rhoc = pcsaft_rhoc(m,σ)\n\nReturns the critical density of the PCSAFT equation of state.\n\nInputs:\n\nm: Segment length (no units)\nσ: Monomer diameter [m]\n\nOutputs:\n\nrhoc : Critical density [mol/m^3]. Returns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_vc","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_vc","text":"vc = pcsaft_vc(m,σ)\n\nReturns the critical volume of the PCSAFT equation of state. Equal to 1/pcsaft_rhoc(m,ϵ,σ)\n\nInputs:\n\nm: Segment length (no units)\nσ: Monomer diameter [m]\n\nOutputs:\n\nvc : Critical volume [m^3/mol]. Returns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#PC-SAFT-saturation-volumes","page":"PC-SAFT","title":"PC-SAFT saturation volumes","text":"","category":"section"},{"location":"pcsaft/","page":"PC-SAFT","title":"PC-SAFT","text":"EoSSuperancillaries.pcsaft_rhosat\nEoSSuperancillaries.pcsaft_vsat\nEoSSuperancillaries.pcsaft_vlsat\nEoSSuperancillaries.pcsaft_vvsat\nEoSSuperancillaries.pcsaft_theta\nEoSSuperancillaries.pcsaft_rhosat_reduced\nEoSSuperancillaries.pcsaft_rholsat_reduced\nEoSSuperancillaries.pcsaft_rhovsat_reduced","category":"page"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_rhosat","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_rhosat","text":"rhol,rhov = pcsaft_rhosat(T,m,ϵ,σ)\n\nReturns the saturation densities of the PCSAFT equation of state at the input temperature T.\n\nInputs:\n\nT: Saturation temperature (Kelvin)\nm: Segment length (no units)\nϵ: Reduced interaction energy [K]\nσ: Monomer diameter [m]\n\nOutputs:\n\nrhol : saturation liquid density [mol/m^3].\nrhov : saturation vapour density [mol/m^3].\n\nReturns NaN,NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_vsat","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_vsat","text":"vl,vv = pcsaft_vsat(T,m,ϵ,σ)\n\nReturns the saturation volumes of the PCSAFT equation of state at the input temperature T.\n\nInputs:\n\nT: Saturation temperature (Kelvin)\nm: Segment length (no units)\nϵ: Reduced interaction energy [K]\nσ: Monomer diameter [m]\n\nOutputs:\n\nvl : saturation liquid volume [m^3/mol].\nvv : saturation vapour volume [m^3/mol].\n\nReturns NaN,NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_vlsat","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_vlsat","text":"vv = pcsaft_vvsat(T,m,ϵ,σ)\n\nReturns the saturation vapour volume of the PCSAFT equation of state at the input temperature T.\n\nInputs:\n\nT: Saturation temperature (Kelvin)\nm: Segment length (no units)\nϵ: Reduced interaction energy [K]\nσ: Monomer diameter [m]\n\nOutputs:\n\nvl : saturation liquid volume [m^3/mol].\n\nReturns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_vvsat","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_vvsat","text":"vl = pcsaft_vlsat(T,m,ϵ,σ)\n\nReturns the saturation liquid volume of the PCSAFT equation of state at the input temperature T.\n\nInputs:\n\nT: Saturation temperature (Kelvin)\nm: Segment length (no units)\nϵ: Reduced interaction energy [K]\nσ: Monomer diameter [m]\n\nOutputs:\n\nvl : saturation liquid volume [m^3/mol].\n\nReturns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and Tmin < T < Tc).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_theta","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_theta","text":"Θ,status = pcsaft_theta(T̃,m)\nΘ,status = pcsaft_theta(T̃,m,T̃c)\n\nCalculates reduced scaled temperature parameter used in PCSAFT saturation volume superancillary\n\nInputs:\n\nT̃: Temperature divided by reduced interaction energy (T/ϵ) (no units)\nm: Segment length (no units)\nT̃c: (Optional) Reduced critical temperature (Tc/ϵ) (no units)\n\nOutputs:\n\nΘ : scaled temperature parameter (no units)\nstatus : Symbol used to signal if the combination of T̃,m is valid or not. It can return one of the following:\n:inrange : if the combination is valid\n:nonfinite : if any input is not finite\n:below_mmin : if m < 1.0\n:over_mmax : if m > 64.0\n:below_Tmin : if T̃ < T̃min, where T̃min = T̃c*exp(-2.20078778)*m^0.37627892\n:over_Tmax : if T̃ > T̃c\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_rhosat_reduced","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_rhosat_reduced","text":"ρ̃l,ρ̃v = pcsaft_rhosat_reduced(Θ,m)\n\nCalculates the reduced saturation densities. The actual densities are calculated as ρᵢ = ρ̃ᵢ/(N_A*σ^3)\n\nInputs:\n\nΘ: reduced temperature parameter (no units)\nm: Segment length (no units)\n\nOutputs:\n\nρ̃l : reduced saturation liquid density (no units)\nρ̃v : reduced saturation vapour density (no units)\n\nReturns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_rholsat_reduced","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_rholsat_reduced","text":"ρ̃l = pcsaft_rholsat_reduced(Θ,m)\n\nCalculates the reduced liquid saturation density. The actual density is ρ = ρ̃/(N_A*σ^3)\n\nInputs:\n\nΘ: reduced temperature parameter (no units)\nm: Segment length (no units)\n\nOutputs:\n\nρ̃l : reduced saturation liquid density (no units)\n\nReturns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).\n\n\n\n\n\n","category":"function"},{"location":"pcsaft/#EoSSuperancillaries.pcsaft_rhovsat_reduced","page":"PC-SAFT","title":"EoSSuperancillaries.pcsaft_rhovsat_reduced","text":"ρ̃v = pcsaft_rhovsat_reduced(Θ,m)\n\nCalculates the reduced vapour saturation density. The actual density is ρ = ρ̃/(N_A*σ^3)\n\nInputs:\n\nΘ: reduced temperature parameter (no units)\nm: Segment length (no units)\n\nOutputs:\n\nρ̃v : reduced saturation vapour density (no units)\n\nReturns NaN if the value is outside the range of the ancillary (1 ≤ m ≤ 64 and 0 < Θ < 1).\n\n\n\n\n\n","category":"function"},{"location":"#EoSSuperancillaries","page":"Home","title":"EoSSuperancillaries","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Superancillary equations to calculate saturation pressures and critical points for selected equations of state.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a Julia port of usnistgov/SAFTsuperanc","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Saturation calculation always involve iterative calculations. those calculations waste a lot of computing power. superancillaries provide direct equations that return some properties within Float64 accuracy.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It can calculate the following properties:","category":"page"},{"location":"","page":"Home","title":"Home","text":"PC-SAFT: critical temperature, critical density, saturated densities.\nvan der Wals: saturated densities, saturated pressures.\nRedlich-Kwong: saturated densities, saturated pressures.\nPeng-Robinson: saturated densities, saturated pressures.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"You can install EoSSuperancillaries by running the following commands on the julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\n\njulia> Pkg.add(\"EoSSuperancillaries\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"you can then load the package into your enviroment by writing in the same REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using EoSSuperancillaries","category":"page"},{"location":"#Cubics","page":"Home","title":"Cubics","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The most known cubic equations of state are defined by two parameters: a and b. The calculation of those parameters depends on each equation of state, values of the fluid at the critical point, and in the case of a, additional parameters that are used to define a temperature dependence. For example, if we take the (canonical) peng robinson equation of state:","category":"page"},{"location":"","page":"Home","title":"Home","text":"a =  Ωa * (R * T)^2 / Pc * α(T)\nb = Ωb * R * T / Pc\nΩa = 0.45723552892138218938\nΩb = 0.077796073903888455972\nα = (1 + κ*(1 - √(T/Tc)))^2\nκ = 0.37464 + 1.54226ω - 0.26992ω²","category":"page"},{"location":"","page":"Home","title":"Home","text":"Where Tc,Pc is the critical temperature and pressure, and ω is the acentric factor. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nThe choice of α(T) does not matter for cubic superancillaries, as we only use the a and b values. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"To use cubic superancillaries, the procedure is the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Calculate a and b values\npass those values to the adequate EoSSuperancillaries method.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using EoSSuperancillaries, Clapeyron\nmodel = PR(\"water\") #model for comparison.\n\n#usage for cubics\na,b,_ = Clapeyron.cubic_ab(model,0.0,373.15)\n#direct calculation using superancillaries, takes about 150 ns\np_eq = pr_psat(373.15,a,b) #96099.13581050883 \nvl,vv = pr_vsat(373.15,a,b)\n#iterative calculation, takes about 2 μs\np_eq2,vl2,vv2 = Clapeyron.saturation_pressure(model,373.15)[1] #96099.13581050835","category":"page"},{"location":"","page":"Home","title":"Home","text":"All cubic superancillaries are implicitly defined in terms of T̃ =RbT/a, and have a range of application of 0.1T̃c < T̃ < T̃c.","category":"page"},{"location":"#PC-SAFT","page":"Home","title":"PC-SAFT","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The parturbed-chain statistical associating fluid theory (PC-SAFT) equation of state is defined, for non-associating compounds, in terms of the segment length (m), the monomer segment diameter (σ) and the reduced interaction energy (ϵ). Associating compounds require additional parameters (associating interaction energy and associating diameter, among with the number and type of sites), and are not supported by the current superancillary present in the package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"One PCSAFT superancillary calculation for the saturation volumes, involves the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"calculate the critical temperature (we can do that with pcsaft_tc(m,ϵ)). in particular T̃c = Tc/ϵ = pcsaft_tc(m,1.0) is used. \ncalculate the minimum temperature (T̃min = T̃c*exp(-2.20078778)*m^0.37627892)\nCalculate Θ = (T̃ - T̃min)/(T̃c - T̃min), where T̃ = T/ϵ\nwith Θ and m, calculate the reduced densities. ρ̃l,ρ̃v\nreescale the densities: ρ = ρ̃/(N_A*σ*σ*σ) where N_A is the Avogadro's constant, equal to 6.02214076e23","category":"page"},{"location":"","page":"Home","title":"Home","text":"All those steps are done automatically by calling pcsaft_vsat (or pcsaft_rhosat), but some steps can be cached. For each compound, you can cache T̃c and T̃min, So you can do the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"function my_pcsaft_rholsat(T,m,ϵ,σ)\n    T̃ = T/ϵ\n    Tc = pcsaft_tc(m,ϵ) #or from cached location\n    T̃c = Tc/ϵ\n    T̃min = T̃c*exp(-2.20078778)*m^0.37627892\n    θ,status = pcsaft_theta(T̃,m,T̃c,T̃min)\n    #you can also check isnan(θ) if you dont care about the exact out of range condition\n    if status == :inrange\n        ρ̃l = pcsaft_rholsat_reduced(θ,m)\n        N_A = 6.02214076e23 #(or EoSSuperancillaries.N_A)\n        return ρ̃l/(N_A*σ*σ*σ)\n    end","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is useful in the cases where you want to benchmark your own saturation solver, or if you are working with extended precision arithmetic.","category":"page"},{"location":"#Notes-on-low-reduced-temperatures","page":"Home","title":"Notes on low reduced temperatures","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"With both cubics and PC-SAFT, the saturation curves are defined in a range Tmin < T < Tc. on low reduced temperatures, the magnitude of the proportion between volumes (vv/vl) reaches the limits of Float64 arithmetic. for solving at temperatures below those ranges, it is recommended to use extended precision numbers (like BigFloat).","category":"page"},{"location":"#References","page":"Home","title":"References","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Bell, I. H., & Deiters, U. K. (2021). Superancillary equations for cubic equations of state. Industrial & Engineering Chemistry Research, 60(27), 9983–9991. doi:10.1021/acs.iecr.1c00847\nBell, I. H., & Deiters, U. K. (2023). Superancillary equations for nonpolar pure fluids modeled with the PC-SAFT equation of state. Industrial & Engineering Chemistry Research. doi:10.1021/acs.iecr.2c02916","category":"page"},{"location":"vdw/","page":"van der Wals","title":"van der Wals","text":"CurrentModule = EoSSuperancillaries","category":"page"},{"location":"vdw/#Contents","page":"van der Wals","title":"Contents","text":"","category":"section"},{"location":"vdw/","page":"van der Wals","title":"van der Wals","text":"Pages = [\"vdw.md\"]","category":"page"},{"location":"vdw/#Index","page":"van der Wals","title":"Index","text":"","category":"section"},{"location":"vdw/","page":"van der Wals","title":"van der Wals","text":"Pages = [\"vdw.md\"]","category":"page"},{"location":"vdw/#van-der-Wals-superancillaries","page":"van der Wals","title":"van der Wals superancillaries","text":"","category":"section"},{"location":"vdw/","page":"van der Wals","title":"van der Wals","text":"EoSSuperancillaries.vdw_vsat\nEoSSuperancillaries.vdw_vlsat\nEoSSuperancillaries.vdw_vvsat\nEoSSuperancillaries.vdw_psat","category":"page"},{"location":"vdw/#EoSSuperancillaries.vdw_vsat","page":"van der Wals","title":"EoSSuperancillaries.vdw_vsat","text":"vl,vv = vdw_vsat(T,a,b)\n\nReturns the saturation liquid and vapour volumes of the vdW cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\nvv : Saturation Vapour Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"vdw/#EoSSuperancillaries.vdw_vlsat","page":"van der Wals","title":"EoSSuperancillaries.vdw_vlsat","text":"vl = vdw_vlsat(T,a,b)\n\nReturns the saturation liquid volume of the vdW cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"vdw/#EoSSuperancillaries.vdw_vvsat","page":"van der Wals","title":"EoSSuperancillaries.vdw_vvsat","text":"vv = vdw_vvsat(T,a,b)\n\nReturns the saturation vapour volume of the vdW cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvv : Saturation Vapour Volume [m^3]. Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"vdw/#EoSSuperancillaries.vdw_psat","page":"van der Wals","title":"EoSSuperancillaries.vdw_psat","text":"p = vdw_psat(T,a,b)\n\nReturns the saturation pressure of the vdW cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\np : Saturation pressure [Pa].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"pr/","page":"Peng-Robinson","title":"Peng-Robinson","text":"CurrentModule = EoSSuperancillaries","category":"page"},{"location":"pr/#Contents","page":"Peng-Robinson","title":"Contents","text":"","category":"section"},{"location":"pr/","page":"Peng-Robinson","title":"Peng-Robinson","text":"Pages = [\"pr.md\"]","category":"page"},{"location":"pr/#Index","page":"Peng-Robinson","title":"Index","text":"","category":"section"},{"location":"pr/","page":"Peng-Robinson","title":"Peng-Robinson","text":"Pages = [\"pr.md\"]","category":"page"},{"location":"pr/#Peng-Robinson-superancillaries","page":"Peng-Robinson","title":"Peng-Robinson superancillaries","text":"","category":"section"},{"location":"pr/","page":"Peng-Robinson","title":"Peng-Robinson","text":"EoSSuperancillaries.pr_vsat\nEoSSuperancillaries.pr_vlsat\nEoSSuperancillaries.pr_vvsat\nEoSSuperancillaries.pr_psat","category":"page"},{"location":"pr/#EoSSuperancillaries.pr_vsat","page":"Peng-Robinson","title":"EoSSuperancillaries.pr_vsat","text":"vl,vv = pr_vsat(T,a,b)\n\nReturns the saturation liquid and vapour volumes of the Peng-Robinson cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\nvv : Saturation Vapour Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"pr/#EoSSuperancillaries.pr_vlsat","page":"Peng-Robinson","title":"EoSSuperancillaries.pr_vlsat","text":"vl = pr_vlsat(T,a,b)\n\nReturns the saturation liquid volume of the Peng-Robinson cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvl : Saturation Liquid Volume [m^3].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"pr/#EoSSuperancillaries.pr_vvsat","page":"Peng-Robinson","title":"EoSSuperancillaries.pr_vvsat","text":"vv = pr_vvsat(T,a,b)\n\nReturns the saturation vapour volume of the Peng-Robinson cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\nvv : Saturation Vapour Volume [m^3]. Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"},{"location":"pr/#EoSSuperancillaries.pr_psat","page":"Peng-Robinson","title":"EoSSuperancillaries.pr_psat","text":"p = pr_psat(T,a,b)\n\nReturns the saturation pressure of the Peng-Robinson cubic equation of state.\n\nInputs:\n\nT: Temperature (Kelvin)\na: Atraction parameter [m^6*Pa/mol]\nb: Covolume [m^3/mol]\n\nOutputs:\n\np : Saturation pressure [Pa].  Returns NaN if the value is outside the range of the ancillary.\n\n\n\n\n\n","category":"function"}]
}
