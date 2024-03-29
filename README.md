# MittagLeffler

[*Mittag-Leffler function*](https://en.wikipedia.org/wiki/Mittag-Leffler_function),

[![Build Status](https://github.com/JuliaMath/MittagLeffler.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaMath/MittagLeffler.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/JuliaMath/MittagLeffler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaMath/MittagLeffler.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![JET QA](https://img.shields.io/badge/JET.jl-%E2%9C%88%EF%B8%8F-%23aa4444)](https://github.com/aviatesk/JET.jl)

```julia
mittlefferr(α,β,z,ρ)   # evaluate Mittag-Leffler function with tolerance ρ
mittlefferr(α,z,ρ)     # mittlefferr(α,1,z,ρ)

mittleff(α,β,z)   # evaluate Mittag-Leffler function with tolerance eps()
mittleff(α,z)     # mittleff(α,1,z)
```

Arguments must satisfy `α > 0`, `β` real, `z` real or complex, `ρ>0`.

For `α<1` and/or `abs(z)<1`, accurate, series-only method are used. The series-only methods work
with BigFloat precision for corresponding input types. Some other parameter ranges also use series
or asymptotic methods.

For some arguments, integrals are evaluated with `quadgk`, with no control on errors. Some results
are accurate, others are not.

```julia
mittleffderiv(α,β,z)   # evaluate derivative of Mittag-Leffler function
mittleffderiv(α,z)     # mittleffderiv(α,1,z)
```

### Bugs

`mittleff` fails for some arguments. In particular, some of those that evaluate integrals.

### Reference

Rudolfo Gorenflo, Joulia Loutchko and Yuri Loutchko,
*Computation of the Mittag-Leffler function and its derivative*,  Fract. Calc. Appl. Anal, **(2002)**
