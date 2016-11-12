using MittagLeffler
using Base.Test

myapp(x,y) = abs(x-y) < 1e-10

@test myapp(mittleff(.5,.5,.5),1.5403698281390346)
@test myapp(mittleff(1.5,.5,.5),1.1448466286155243)
@test myapp(mittleff(2.3, .7 + 2. * im), 1.201890136368392 + 0.7895394560075035 * im)
@test myapp(mittleff(2.3, .7 + 0.2 * im) , (1.268233154873853 + 0.07914994421659409im))
nothing
