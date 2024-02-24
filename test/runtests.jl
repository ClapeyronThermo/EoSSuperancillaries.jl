using EoSSuperancillaries
using Test
const ES = EoSSuperancillaries
using EoSSuperancillaries: cheb_eval

@testset "cubic" begin
    T̃ = 0.125
    @testset "vdW" begin
        @test cheb_eval(ES._vdw_p(),T̃) ≈ 0.0002958543239347111
        @test cheb_eval(ES._vdw_vl(),T̃) ≈ 0.8536251284168529
        @test cheb_eval(ES._vdw_vv(),T̃) ≈ 0.002407389267319304
    end

    @testset "RK" begin
        @test cheb_eval(ES._rk_p(),T̃) ≈ 0.001736846506201768
        @test cheb_eval(ES._rk_vl(),T̃) ≈ 0.6976615743280177
        @test cheb_eval(ES._rk_vv(),T̃) ≈ 0.01555500889873714
    end

    @testset "PR" begin
        @test cheb_eval(ES._pr_p(),T̃) ≈ 0.003034198868923775
        @test cheb_eval(ES._pr_vl(),T̃) ≈ 0.6394564580846998
        @test cheb_eval(ES._pr_vv(),T̃) ≈ 0.03023195086998487
    end
    # Write your tests here.
end
#=
1 1.1481738529594 0.578301305568106 0.081366617928304
1 0.893024107857271 0.742270101982426 0.014085538623874
64 3.79621938847095 0.00606038110809739 5.0222327626229e-06
64 2.952615079922 0.00927005575174636 3.25402208337432e-11

=#
@testset "PCSAFT" begin
    m1,m64 = 1,64
    T̃c1 = ES.pcsaft_tc(m1,1.0)
    T̃c64 = ES.pcsaft_tc(m64,1.0)
    theta_09_01 = ES.pcsaft_theta(0.9*T̃c1,m1)[1]
    theta_07_01 = ES.pcsaft_theta(0.7*T̃c1,m1)[1]
    theta_09_64 = ES.pcsaft_theta(0.9*T̃c64,m64)[1]
    theta_07_64 = ES.pcsaft_theta(0.7*T̃c64,m64)[1]
    @test 0.9*T̃c1 ≈ 1.1481738529594
    rhol_m1_09,rhov_m1_09 = ES.pcsaft_rhosat_reduced(theta_09_01,m1)
    rhol_m1_07,rhov_m1_07 = ES.pcsaft_rhosat_reduced(theta_07_01,m1)
    rhol_m64_09,rhov_m64_09 = ES.pcsaft_rhosat_reduced(theta_09_64,m64)
    rhol_m64_07,rhov_m64_07 = ES.pcsaft_rhosat_reduced(theta_07_64,m64)
    @test rhol_m1_09 ≈ 0.578301305568106
    @test rhov_m1_09 ≈ 0.081366617928304
    @test rhol_m1_07 ≈ 0.742270101982426
    @test rhov_m1_07 ≈ 0.014085538623874
    @test rhol_m64_09 ≈ 0.00606038110809739
    @test rhov_m64_09 ≈ 5.0222327626229e-06
    @test rhol_m64_07 ≈ 0.00927005575174636
    @test rhov_m64_07 ≈ 3.25402208337432e-11
    @test pcsaft_vc(1.0,3.7039e-10) ≈ 0.00010836057852600164
end