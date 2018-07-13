## This code implements algorithms in
## Rudolfo Gorenflo, Joulia Loutchko and Yuri Loutchko,
## *Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal, **(2002)**

# When support for v0.6 is dropped, use fully qualified name
if VERSION >= v"0.7-"
    import SpecialFunctions: gamma
end

macro br(n)
#    esc(:(println(STDERR, "branch ", $n)))
    nothing
end

function P(α,β,ϵ,ϕ,z)
    ω = ϕ * (1+(1-β)/α) + ϵ^(1/α) * sin(ϕ/α)
    return (1/(2*α*pi)) * ϵ^(1+(1-β)/α)*exp(ϵ^(1/α) * cos(ϕ/α)) * (cos(ω) + im * sin(ω))/(ϵ*exp(im*ϕ)-z)
end

if VERSION >= v"0.5-"
    import QuadGK: quadgk
end

ourquadgk(f,a,b) = quadgk(f,a,b; order=7)[1]

function Kint(α,β,a,χ0,z)
    return ourquadgk(χ -> K(α,β,χ,z), a, χ0)
end

function K(α,β,χ,z)
    den = (χ^2-2*χ*z*cos(α*pi)+z^2)
    return (1/(α*pi)) * χ^((1-β)/α) * exp(-χ^(1/α))*( χ * sin(pi*(1-β)) - z * sin(pi*(1-β+α)))/den
end

Kint(α,β,χ0,z) = Kint(α,β,0,χ0,z)

mpow(x::Complex, y) = x^y
mpow(x::Real, y) = x >= 0 ? x^y : Complex(x,0)^y

# The following gymnastics are to get around the type-unstable mpow
function sum2(α, β, z::Real, k0)
    return z > 0 ? sum2_pos(α, β, z, k0) : sum2_neg(α, β, z, k0)
end
sum2(α, β, z::Complex, k0) = sum2_neg(α, β, z, k0)

for funcname in (:sum2_neg, :sum2_pos)
    local zarg
    local stype
    if funcname == :sum2_neg
        zarg = :( (z + Complex(0,0) ))
        stype = :( typeof(z + Complex(0,0) ) )
    else
        zarg = :z
        stype = :( typeof(z) )
    end
    @eval function ($funcname)(α, β, z, k0)
        s::($stype) = zero($stype)
        for k=1:k0
            arg = β - α * k
            if !( round(arg) == arg && arg < 0)
                s += ($zarg)^(-k) / gamma(arg)
            end
        end
        return s
    end
end

function choosesum(α,β,z,ρ)
    k0 = floor(Int, -log(ρ)/log(abs(z)))
    if abs(angle(z)) < pi*α/4 + 1//2 * min(pi,pi*α)
        @br 3
        return 1/α * z^((1-β)/α) * exp(z^(1/α)) - sum2(α,β,z,k0)
    else
        @br 4
        return - sum2(α,β,z,k0)
    end
end


function mittleffsum(α, β, z)
    @br 1
    k0::Int = floor(Int, α) + 1
    s = zero(z)
    for k=0:(k0 - 1)
        s += mittleff(α / k0, β, mpow(z, (1 // k0)) * exp(2pi * im * k / k0))
    end
    return s / k0
end

function mittleffsum2(α,β,z,ρ)
    @br 2
    k0 = max(ceil(Int,(1-β)/α), ceil(Int, log(ρ*(1-abs(z)))/log(abs(z))))
    s = zero(z)
    for k=0:k0
        s += z^k/gamma(β+α*k)
    end
    return s
end

Pint(α,β,ϵ,z) = ourquadgk(ϕ -> P(α,β,ϵ,ϕ,z), -α*pi, α*pi)
Pint(α,β,z) = Pint(α,β,1,z)

function mittleffints(α,β,z,ρ)
    az = abs(z)
    ab = abs(β)
    χ0 = β >= 0 ?
      max(1,2*az,(-log(pi*ρ/6))^α) :
      max((ab+1)^α, 2*az,(-2*log( pi*ρ/(6*(ab+2)*(2*ab)^ab)))^α)
    aaz = abs(angle(z))
    if aaz > α * pi
        if β <= 1
            @br 5
            return Kint(α,β,χ0,z)
        else
            @br 6
            return Kint(α,β,1,χ0,z) + Pint(α,β,z)
        end
    elseif aaz < pi*α
        if β <= 1
            @br 7
            return Kint(α,β,χ0,z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        else
            @br 8
            return Kint(α,β,az/2,χ0,z) + Pint(α,β,(az)/2,z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        end
    else
        @br 9
        return Kint(α,β,(az+1)/2,χ0,z) + Pint(α,β,(az+1)/2,z)
    end
end

"""
    mittlefferr(α,z,ρ)

Compute mittlefferr(α,1,z,ρ).
"""
mittlefferr(α,z,ρ) = mittlefferr(α,1,z,ρ)


"""
    mittlefferr(α,β,z,ρ)

Compute the Mittag-Leffler function at `z` for parameters `α,β` with
accuracy `ρ`.
"""
function mittlefferr(α,β,z,ρ)
    ρ > 0 || throw(DomainError())
    _mittlefferr(α,β,z,ρ)
end

_mittlefferr(α::Real,β::Real,z::Real,ρ::Real) = real(_mittleff(α,β,z,ρ))
_mittlefferr(α::Real,β::Real,z::Complex,ρ::Real) = _mittleff(α,β,z,ρ)

# The second definition would work for both complex and real
myeps(x) = x |> one |> float |> eps
myeps(x::Complex) =  x |> real |> myeps


"""
    mittleff(α,β,z)

Compute the Mittag-Leffler function at `z` for parameters `α,β`.
"""
mittleff(α, β, z) = _mittlefferr(α,β,z,myeps(z))
mittleff(α, β, z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, β, float(z))

"""
    mittleff(α,z)

Compute `mittleff(α,1,z)`.
"""
mittleff(α,z) = _mittlefferr(α,1,z,myeps(z))
mittleff(α,z::Union{Integer,Complex{T}}) where {T<:Integer} = mittleff(α, float(z))

function _mittleff(α,β,z,ρ)
#    if β == 1
        # if α == 1/2   # disable this. There is never an error here. But, this triggers a mysterious, untraceable bug in quadgk
        #     res = try
        #         exp(z^2)*erfc(-z)
        #     catch
        #         error("Failed in exp. erfc")
        #     end
        #     return res
        # end
        # Disable these, because we would have to take care of domain errors, etc.
        # α == 0 && return 1/(1-z)
        # α == 1 && return exp(z)
        # α == 2 && return cosh(sqrt(z))
        # α == 3 && return (1//3)*(exp(z^(1//3)) + 2*exp(-z^(1//3)/2) * cos(sqrt(convert(typeof(z),3))/2 * z^(1//3)))
        # α == 4 && return (1//2)*(cosh(z^(1//4)) + cos(z^(1//4)))
#    end
    z == 0 && return 1/gamma(β)
    α == 1 && β == 1 && return(exp(z))
    α < 0  && throw(DomainError())
    az = abs(z)
    1 < α && return mittleffsum(α,β,z)
    az < 1 && return mittleffsum2(α,β,z,ρ)
    az > floor(10+5*α) && return choosesum(α,β,z,ρ)
    return mittleffints(α,β,z,ρ)
end

_mittleff(α,β,z) = mittleff(α,β,z,myeps(z))


"""
    mittleffderiv(α,β,z)

Compute the derivative of the Mittag-Leffler function at `z` for parameters `α,β`.
"""
function mittleffderiv(α, β, z)
    #derivative of Mittag Leffler function WRT to main argument Z.
    #take q = 0.5. Paper requires |z| <= q < 1.
    q = 1//2
    #case 1, small z
    if abs(z) <= q
        ω = α + β - 3//2
        D = α^2 - 4*α*β + 6*α + 1
        #k1
        if α>1
            k₁ = ((2-α-β)/(α-1)) + 1
        elseif α>0 && α<=1 && D<=0
            k₁ = ((3-α-β)/α) + 1
        else
            k₁ = maximum([((3-α-β)/α) + 1, ((1-2*ω*α+sqrt(D))/(2*(α^2)))+1])
        end
        k₀ = maximum([k₁, log(myeps(z)*(1-abs(z)))/log(abs(z))])
        k₀ = ceil(Int,k₀) #take ceiling (not specified in paper whether floor or ceiling)
        out = zero(z)
        for k in 0:k₀
            out = out + ((k+1)*z^k)/(gamma(α+β+α*k))
        end
    #case 2, larger z
    else
        out = (mittleff(α,β-1,z) - (β-1)*mittleff(α,β,z))/(α*z)
    end
    return out
end

"""
    mittleffderiv(α,z)

Compute mittleffderiv(α,1,z)
"""
mittleffderiv(α, z) = mittleffderiv(α,1,z)
