
function P(α,β,ϵ,ϕ,z)
    ω = ϕ * (1+(1-β)/α) + ϵ^(1/α) * sin(ϕ/α)    
    (1/(2*α*pi)) * ϵ^(1+(1-β)/α)*exp(ϵ^(1/α) * cos(ϕ/α)) * (cos(ω) + im * sin(ω))/(ϵ*exp(im*ϕ)-z)
end

# TODO: keep track of error
ourquadgk(f,a,b) = quadgk(f,a,b)[1]

Pint(α,β,z) = Pint(α,β,1,z)
Pint(α,β,ϵ,z) = ourquadgk( ϕ -> P(α,β,1,ϕ,z), -α*pi, α*pi)

K(α,β,χ,z) = (1/(α*pi)) * χ^((1-β)/α) * exp(-χ^(1/α))*( χ * sin(pi*(1-β)) - z * sin(pi*(1-β-α)))/(χ^2-2*χ*z*cos(α*pi)+z^2)

Kint(α,β,χ0,z) = Kint(α,β,0,χ0,z)

Kint(α,β,a,χ0,z) = ourquadgk(χ -> K(α,β,χ,z), a, χ0)

function mittleffsum(α,β,z)
    k0 = floor(Int,α) + 1
    s = zero(z)
    for k=0:(k0-1)
        s += mittleff(α/k0,β,z^(1/k0)*exp(2*pi*im*k/k0))
    end
    s / k0
end

function mittleffsum2(α,β,z,ρ)
    k0 = max(ceil(Int,(1-β)/α), ceil(Int, log(ρ*(1-abs(z)))/log(abs(z))))
    s = zero(z)
    for k=0:k0
        s += z^k/gamma(β+α*k)
    end
    s
end

function sum2(α,β,z,k0)
    s = zero(z)
    for k=1:k0
        s += z^(-k)/gamma(β-α * k)
    end
    s
end

function choosesum(α,β,z,ρ)
    k0 = floor(Int, -log(ρ)/log(abs(z)))
    if abs(angle(z)) < pi*α/4 + 1//2 * min(pi,pi*α)
        return 1/α * z^((1-β)/α) * exp(z^(1/α)) - sum2(α,β,z,k0)
    else
        return - sum2(α,β,z,k0)
    end
end

function mittleffints(α,β,z,ρ)
    az = abs(z)
    ab = abs(β)
    χ0 = β >= 0 ?
      max(1,2*az,(-log(pi*ρ/6))^α) :
      max((ab+1)^α, 2*az,(-2*log( pi*ρ/(6*(ab+2)*(2*ab)^ab)))^α)
    aaz = abs(angle(z))
    if aaz > α * pi
        if β <= 1
            return Kint(α,β,χ0,z)
        else
            return Kint(α,β,1,χ0,z) + Pint(α,β,z)
        end
    elseif aaz < pi*α
        if β <= 1
            return Kint(α,β,χ0,z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        else
            return Kint(α,β,abs(z)/2,χ0,z) + Pint(α,β,(abs(z))/2,z) + (1/α)*z^((1-β)/α) * exp(z^(1/α))
        end
    else
        return Kint(α,β,(abs(z)+1)/2,χ0,z) + Pint(α,β,(abs(z)+1)/2,z)
    end
end

## TODO: Do real values sometimes return complex result ?
mittlefferr(α,z,ρ) = mittlefferr(α,1,z,ρ)
mittlefferr(α::Real,β::Real,z::Real,ρ::Real) = real(_mittleff(α,β,z,ρ))
mittlefferr(α::Real,β::Real,z::Complex,ρ::Real) = _mittleff(α,β,z,ρ)

mittleff(α,β,z) = mittlefferr(α,β,z,eps())
mittleff(α,z) = mittlefferr(α,1,z,eps())

function _mittleff(α,β,z,ρ)
    1 < α && return mittleffsum(α,β,z)
    z == 0 && return 1/gamma(β)
    abs(z) < 1 && return mittleffsum2(α,β,z,ρ)
    abs(z) > floor(10+5*α) && return choosesum(α,β,z,ρ)
    mittleffints(α,β,z,ρ)
end

_mittleff(α,β,z) = mittleff(α,β,z,eps())
