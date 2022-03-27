function mittleff_matlab(α,β,γ,z)
    (α,β,γ,z)=promote(α,β,γ,z)
    T=typeof(α)
    if real(α)<=0 || real(γ)<=0 || !isreal(α) || !isreal(β) || !isreal(γ)
        ArgumentError("Error in the parameters of the Mittag-Leffler function. α($α) and γ($γ) must be real and positive. β($β) must be real.")
    end
    if abs(γ-1) > eps(T)
        if α>1
            ArgumentError("With the three parameters Mittag-Leffler function, the parameter α($α) must satisfy 0<α<1.")
        end
        if abs(angle(z*(abs(z)>eps(T)))) <= α*π
            ArgumentError("With the three parameters Mittag-Leffler function, this code works only when |arg(z)|>απ. z=$z")
        end
    end
    ϵ=10*eps(T)
    if abs(z)<ϵ
        return 1/gamma(β)
    else
        return LTInversion(one(z),z,α,β,γ,ϵ)
    end
end
mittleff_matlab(α,β,z)=mittleff_matlab(α,β,1,z)
mittleff_matlab(α,z)=mittleff_matlab(α,1,z);

function LTInversion(t,λ,α,β,γ,ϵ)
    T=typeof(λ)

    # Evaluation of the relevant poles
    θ = angle(λ)
    kmin = ceil(-α/2-θ/2/π)
    kmax = floor(α/2-θ/2/π)
    k_vett = kmin:kmax
    s_star = abs(λ)^(1/α)*exp.(1im*(θ.+2π*k_vett)./α)

    # Evaluation of ϕ(s_star) for each pole
    ϕ_s_star = (real.(s_star)+abs.(s_star))/2

    # Sorting of the poles according to the value of ϕ(s_star)
    index_s_star = sortperm(ϕ_s_star)
    ϕ_s_star = ϕ_s_star[index_s_star]
    s_star = s_star[index_s_star]

    # Deleting possible poles with ϕ_s_star=0
    index_save = ϕ_s_star.>ϵ
    s_star = s_star[index_save]
    ϕ_s_star = ϕ_s_star[index_save]

    # Inserting the origin in the set of the singularities
    s_star = vcat(0,s_star)
    ϕ_s_star = vcat(0,ϕ_s_star)
    J1 = length(s_star); J=J1-1

    # Strength of the singularities
    p = vcat(max(0,-2*(α*γ-β+1)), fill(γ,(J,)))
    q = vcat(fill(γ,(J,)), T(Inf))
    append!(ϕ_s_star, T(Inf))

    # Looking for the admissible regions with respect to round-off errors
    admissible_regions = findall((ϕ_s_star[1:end-1].<(log(ϵ)-log(eps(T)))/t) .& (ϕ_s_star[1:end-1].<ϕ_s_star[2:end]))

    # Initializing vectors for optimal parameters
    JJ1 = admissible_regions[end]
    μ_vett = fill(T(Inf),JJ1)
    N_vett = fill(T(Inf),JJ1)
    h_vett = fill(T(Inf),JJ1)

    # Evaluation of parameters for inversion of LT in each admissible region
    find_region=false
    while !find_region
        for j1 = admissible_regions
            if j1<J1
                (μj,hj,Nj) = OptimalParam_RB(t,ϕ_s_star[j1],ϕ_s_star[j1+1],p[j1],q[j1],ϵ)
            else
                (μj,hj,Nj) = OptimalParam_RU(t,ϕ_s_star[j1],p[j1],ϵ)
            end
            μ_vett[j1]=μj; h_vett[j1]=hj; N_vett[j1]=Nj
        end
        if minimum(N_vett)>200
            ϵ *= 10
        else
            find_region = true
        end
    end

    # Selection of the admissible region for integration which involves the minimum number of nodes 
    (N,iN)=findmin(N_vett); μ=μ_vett[iN]; h=h_vett[iN]

    # Evaluation of the inverse Laplace transform
    k=-N:N; u=h*k
    z = μ*(im*u.+1).^2
    zd = -2*μ*u .+ 2im*μ
    zexp = exp.(t*z)
    F = z.^(α*γ-β) ./ (z.^α.-λ).^γ .* zd
    S = zexp .* F;
    Integral = h*sum(S)/(2*im*π);

    # Evaluation of residues
    ss_star = s_star[iN+1:end];
    Residues = sum(1/α * ss_star.^(1-β) .* exp.(t*ss_star));

    # Evaluation of the ML function
    E = Integral + Residues;
    if isreal(λ)
        E = real(E)
    end
    return E;
end

function OptimalParam_RB(t,ϕ_s_star_j,ϕ_s_star_j1,pj,qj,ϵ)
    # Definition of some constants
    T=typeof(t); fac=1.01; conservative_error_analysis=false;

    # Maximum value of fbar as the ration between tolerance and round-off unit
    f_max = ϵ/eps(T);

    # Evaluation of the starting values for sq_ϕ_star_j and sq_ϕ_star_j1
    sq_ϕ_star_j = sqrt(ϕ_s_star_j);
    threshold = 2*sqrt(log(f_max)/t);
    sq_ϕ_star_j1 = min(sqrt(ϕ_s_star_j1), threshold-sq_ϕ_star_j);

    # Zero or negative values of pj and qj
    if pj<10*ϵ && qj<10*ϵ
        sq_φ_star_j=sq_ϕ_star_j; sq_φ_star_j1=sq_ϕ_star_j1; adm_region=true;
    # Zero or negative values of just pj
    elseif pj<10*ϵ && qj>=10*ϵ
        sq_φ_star_j=sq_ϕ_star_j;
        if sq_ϕ_star_j>0
            f_min = fac*(sq_ϕ_star_j/(sq_ϕ_star_j1-sq_ϕ_star_j))^qj
        else
            f_min = fac;
        end
        if f_min<f_max
            f_bar = f_min+f_min/f_max*(f_max-f_min);
            fq = f_bar^(-1/qj);
            sq_φ_star_j1 = (2*sq_ϕ_star_j1-fq*sq_ϕ_star_j)/(2+fq);
            adm_region=true;
        else
            adm_region=false;
        end
    # Zero or negative values of just qj
    elseif pj>=10*ϵ && qj<10*ϵ
        sq_φ_star_j1=sq_ϕ_star_j1;
        f_min = fac*(sq_ϕ_star_j1/(sq_ϕ_star_j1-sq_ϕ_star_j))^qj
        if f_min<f_max
            f_bar = f_min+f_min/f_max*(f_max-f_min);
            fp = f_bar^(-1/pj);
            sq_φ_star_j = (2*sq_ϕ_star_j+fp*sq_ϕ_star_j1)/(2-fp);
            adm_region=true;
        else
            adm_region=false;
        end
    # Positive values of both pj and qj
    elseif pj>=10*ϵ && qj>=10*ϵ
        f_min = fac*(sq_ϕ_star_j+sq_ϕ_star_j1)/(sq_ϕ_star_j1-sq_ϕ_star_j)^max(pj,qj);
        if f_min<f_max
            f_min=max(f_min,1.5)
            f_bar = f_min+f_min/f_max*(f_max-f_min);
            fp = f_bar^(-1/pj);
            fq = f_bar^(-1/qj);
            if !conservative_error_analysis
                w = -ϕ_s_star_j1*t/log(ϵ);
            else
                w = -2*ϕ_s_star_j1*t/(log(ϵ)-ϕ_s_star_j1*t)
            end
            den = 2+w-(1+w)*fp+fq;
            sq_φ_star_j = ((2+w+fq)*sq_ϕ_star_j+fp*sq_ϕ_star_j1)/den;
            sq_φ_star_j1 = (-(1+w)*fq*sq_ϕ_star_j+(2+w-(1+w)*fp)*sq_ϕ_star_j1)/den;
            adm_region=true;
        else
            adm_region=false;
        end
    end

    if adm_region
        ϵ /= f_bar;
        if !conservative_error_analysis
            w = -sq_φ_star_j1^2*t/log(ϵ)
        else
            w = -2*sq_φ_star_j1^2*t/(log(ϵ)-sq_φ_star_j1^2*t)
        end
        μj = (((1+w)*sq_φ_star_j+sq_φ_star_j1)/(2+w))^2;
        hj = -2π/log(ϵ)*(sq_φ_star_j1-sq_φ_star_j)/((1+w)*sq_φ_star_j+sq_φ_star_j1)
        Nj = ceil(sqrt(1-log(ϵ)/t/μj)/hj);
    else
        μj=zero(T); hj=zero(T); Nj=T(Inf);
    end
    return (μj,hj,Nj)
end

function OptimalParam_RU(t,ϕ_s_star_j,pj,ϵ)
    # Evaluation of the starting values for sq_phi_star_j
    sq_ϕ_s_star_j=sqrt(ϕ_s_star_j)
    if ϕ_s_star_j>0
        φ_star_j=ϕ_s_star_j*1.01
    else
        φ_star_j=0.01;
    end
    sq_φ_star_j=sqrt(φ_star_j);

    # Definition of some constants
    f_min=1; f_max=10; f_tar=5;
    T=typeof(t); log_ϵ=log(ϵ);

    # Iterative process to look for fbar in [f_min,f_max]
    stop=false; sq_μj=0; A=0;
    while !stop
        ϕ_t=φ_star_j*t; log_ϵ_ϕ_t=log_ϵ/ϕ_t;
        Nj = ceil(ϕ_t/π*(1-3*log_ϵ_ϕ_t/2+sqrt(1-2*log_ϵ_ϕ_t)))
        A = π*Nj/ϕ_t
        sq_μj = sq_φ_star_j*abs(4-A)/abs(7-sqrt(1+12*A))
        fbar = ((sq_φ_star_j-sq_ϕ_s_star_j)/sq_μj)^(-pj)
        stop = (pj<10*ϵ) || (f_min<fbar && fbar<f_max)
        if !stop
            sq_φ_star_j = f_tar^(-1/pj)*sq_μj+sq_ϕ_s_star_j
            φ_star_j = sq_φ_star_j^2
        end
    end
    μj=sq_μj^2
    hj = (-3*A-2+2*sqrt(1+12*A))/(4-A)/Nj;

    # Adjusting integration parameters to keep round-off errors under control
    log_eps=log(eps(T)); threshold=(log_ϵ-log_eps)/t;
    if μj>threshold
        if abs(pj)<10*ϵ
            Q=0
        else
            Q=f_tar^(-1/pj)*sqrt(μj)
        end
        if φ_star_j<threshold
            w = sqrt(log_eps/(log_eps-log_ϵ))
            u = sqrt(-φ_star_j*t/log_eps)
            μj= threshold
            Nj= ceil(w*log_ϵ/(2π*(u*w-1)))
            hj= w/Nj
        else
            Nj=T(Inf); hj=zero(T);
        end
    end
    return (μj,hj,Nj)
end
