using MittagLeffler
using Aqua: Aqua

const ThePackage = MittagLeffler

if VERSION >= v"1.1"
    @testset "aqua piracies" begin
        Aqua.test_piracies(ThePackage)
    end
end

@testset "aqua deps compat" begin
    Aqua.test_deps_compat(ThePackage)
end

@testset "aqua unbound_args" begin
    Aqua.test_unbound_args(ThePackage)
end

@testset "aqua undefined exports" begin
    Aqua.test_undefined_exports(ThePackage)
end

# Perhaps some of these should be fixed. Some are for combinations of types
# that make no sense.
@testset "aqua test ambiguities" begin
    Aqua.test_ambiguities([ThePackage, Core, Base])
end

@testset "aqua project extras" begin
    Aqua.test_project_extras(ThePackage)
end

@testset "aqua state deps" begin
    Aqua.test_stale_deps(ThePackage)
end
