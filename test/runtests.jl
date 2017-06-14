using MCMCDiagnostics
using Base.Test

@testset "IID samples" begin
    chains = [randn(1000) for _ in 1:10]
    stats = convergence_statistics.(chains)

    @test 9000 ≤ effective_sample_size(stats...) ≤ 10000

    @test 1 ≤ potential_scale_reduction(stats...) ≤ 1.001
end

@testset "AR(1)" begin
    function ar1(ρ, N; σ=1.0)
        x = Vector{Float64}(N)
        z = σ*randn()/√(1-ρ^2)
        for i in 1:N
            z = ρ*z + randn()*σ
            x[i] = z
        end
        x
    end
    "Theoretical ESS."
    ar1_ess_factor(ρ) = 1/(1+2*ρ/(1-ρ))
    ess(ρ, N) = effective_sample_size(ar1(ρ, N, σ=randn()^2+0.5))
    @test
end
