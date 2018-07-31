using MittagLeffler
using Test

@testset "mittleff" begin
    @test isapprox(mittleff(0.5, 0.5, 0.5),  1.5403698281390346)
    @test isapprox(mittleff(1.5, 0.5, 0.5), 1.1448466286155243)
    @test isapprox(mittleff(2.3, 0.7 + 2.0 * im), 1.201890136368392 + 0.7895394560075035 * im)
    @test isapprox(mittleff(2.3, 0.7 + 0.2 * im) , (1.268233154873853 + 0.07914994421659409im))
    @test mittleff(big".3",100.0) ==
        big"8.721285946907744692995882256235296113802695745418015206361825134909144332670706e+2015816"
    # allow alpha == beta == 1  --> exp(z)
    @test isapprox(mittleff(1, -2), 0.1353352832366127)
end

@testset "branches" begin
    # convert z to float or complex float to avoid integer to negative power.
    @test isapprox(mittleff(0.9, 0.5, 22 + 22im), -2.7808021618204008e13 - 2.8561425165239754e13im)
    # Test branch that calls `Pint`
    @test isapprox(mittleff(0.1, 1.05, 0.9 + 0.5im), 0.17617901349590603 + 2.063981943021305im)
end

@testset "derivative" begin
    myder(f,x,h) = (f(x+h/2)-f(x-h/2))/h
    @test isapprox(myder(z -> mittleff(.4,z),.4,1e-5), mittleffderiv(.4,.4))
    @test abs(myder(z -> mittleff(big".4",z),big".4",BigFloat(1//10^16)) - mittleffderiv(big".4",big".4")) < 1e-32
end

