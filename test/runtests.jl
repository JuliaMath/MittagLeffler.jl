using MittagLeffler
using Base.Test

myapp(x,y) = abs(x-y) < 1e-10

@test myapp(mittleff(.5,.5,.5),1.5403698281390346)
@test myapp(mittleff(1.5,.5,.5),1.1448466286155243)
@test myapp(mittleff(2.3, .7 + 2. * im), 1.201890136368392 + 0.7895394560075035 * im)
@test myapp(mittleff(2.3, .7 + 0.2 * im) , (1.268233154873853 + 0.07914994421659409im))
@test mittleff(big".3",100.0) == big"8.721285946907744692995882256235296113802695745418015206361825134909144332670706e+2015816"

# allow alpha == beta == 1  --> exp(z)
@test myapp(mittleff(1,-2), 0.1353352832366127)

myder(f,x,h) = (f(x+h/2)-f(x-h/2))/h

@test myapp(myder(z -> mittleff(.4,z),.4,1e-5), mittleffderiv(.4,.4))
@test abs(myder(z -> mittleff(big".4",z),big".4",BigFloat(1//10^16)) - mittleffderiv(big".4",big".4")) < 1e-32

nothing
