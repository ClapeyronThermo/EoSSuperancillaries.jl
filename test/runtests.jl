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
