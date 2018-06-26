using MCMCDiagnostics

using StatsBase: var

# TODO: change/remove when dropping support for v0.7
using Compat
using Compat.Test
using Compat.Random: rand, randn, srand
range = Compat.range

# consistent testing
srand(UInt32[0x3d50884f, 0xd6560f94, 0x6c04ab37, 0xb1c52878])

@testset "IID samples" begin
    chains = [randn(1000) for _ in 1:10]
    @test 9900 ≤ effective_sample_size(vcat(chains...)) ≤ 10100
    @test 1 ≤ potential_scale_reduction(chains...) ≤ 1.001
end

@testset "IID dispersed ragged PSRF" begin
    μs = -10:1:10
    chains = [randn(1000) * rand(1:10) .+ μ for μ in μs]
    B = var(μs)
    W = mean(var.(chains))
    expected_R̂ = √(1 + B / W)
    R̂ = potential_scale_reduction(chains...)
    @test 0.99 ≤ (R̂ / expected_R̂) ≤ 1.01
end

@testset "AR(1)" begin
    "Simulate an AR(1) process (coefficient ρ, N draws, std σ)."
    function simulate_ar1(ρ, N; σ = 1.0)
        x = Vector{Float64}(undef, N)
        z = σ * randn() / √(1 - ρ^2)
        for i in 1:N
            z = ρ*z + randn()*σ
            x[i] = z
        end
        x
    end

    "Theoretical ESS."
    ar1_ess_factor(ρ) = 1/(1 + 2*ρ/(1-ρ))

    """
    Estimate of effective sample size for an AR(1) process with
    coefficient ρ, divided by the theoretical value.  Should be
    invariant to σ, and around 1.
    """
    function rel(ρ, N)
        effective_sample_size(simulate_ar1(ρ, N, σ=randn()^2+0.5))/(ar1_ess_factor(ρ)*N)
    end

    for ρ in range(0.1; length = 10, stop = 0.4)
        @test 0.95 ≤ rel(ρ,100000) ≤ 1.05
    end

    for ρ in range(0.5; length = 10, stop = 0.9)
        @test 0.95 ≤ rel(ρ,1000000) ≤ 1.05
    end

    # testing cap, since for ρ < 0, ESS ≥ N
    for ρ in range(-0.5; length = 10, stop = -0.9)
        effective_sample_size(simulate_ar1(ρ, 1000, σ=randn()^2+0.5)) == 1000
    end
end
