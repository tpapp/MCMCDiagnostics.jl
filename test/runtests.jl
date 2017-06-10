using MCMCDiagnostics
using Base.Test

@testset "IID samples" begin
    chains = [randn(1000) for _ in 1:10]
    stats = chain_statistics.(chains)

    @test 9000 ≤ effective_sample_size(stats...) ≤ 10000

    @test 1 ≤ potential_scale_reduction(stats...) ≤ 1.001
end
