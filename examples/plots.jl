using MittagLeffler
using Plots

function mittleffplot(α,β,t0,t1)
    plot( x ->  mittleff(α,β,-x), t0,t1)
end

function mittleffplot2(α,β,t0,t1,arg; kws...)
    st = (t1-t0)/20
    plot( r ->  abs(mittlefferr(α,β,r*exp(im*α*arg), 10*eps())), t0:st:t1; kws...)
end

function mittleffplot3(α,β,t0,t1)
    st = (t1-t0)/1000
    plot( r ->  abs(mittleff(α,β,r*exp(im*α*pi/2))), t0:st:t1, ylims = (1.3,1.38))
end
