@testset "Distribution functions" begin

    @test Boltzmann(0.1) ≈ 1.1051709180756477
    @test Blakemore(0.1) ≈ 0.8511816057678475
    @test FermiDiracMinusOne(0.1) ≈ 0.52497918747894
    @test FermiDiracOneHalfBednarczyk(0.1) ≈ 0.8277989807317992
    @test FermiDiracOneHalfTeSCA(0.1) ≈ 0.833056882078161
    @test GaussFermiPaasch(2.0)(0.1) ≈ 0.5148909449738751
    @test GaussFermiSimpson13(2.0, 0.0, 1.0, 10_000)(0.1) ≈ 0.5151388501936096

end
