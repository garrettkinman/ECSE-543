using Test

"""
Uses Choleski decomposition to solve `𝐀𝐱 = 𝐛`, where 𝐀 is real, symmetric, and positive-definite. Returns the vector 𝐱.
"""
function choleski(𝐀::AbstractMatrix{<:Real}, 𝐛::AbstractVector{<:Real}, halfband=nothing)
    n, m = size(𝐀)

    # simple error checking
    if n != m
        throw(DimensionMismatch("Matrix 𝐀 must be square"))
    elseif length(𝐛) != n
        throw(DimensionMismatch("Matrix 𝐀 must be n×n, and vector 𝐛 must be n×1"))
    end
    
    # initialize
    𝐋 = zeros(n,n)
    𝐲 = zeros(n,1)
    𝐱 = zeros(n,1)

    #decompose 𝐋
    for j ∈ 1:n
        sumⱼ = 0
        for q = 1:j-1
            @inbounds sumⱼ = sumⱼ + 𝐋[j,q]^2
        end
        @inbounds 𝐋[j,j] = sqrt(𝐀[j,j] - sumⱼ)

        sumᵢ = 0

        i_range = (j+1):n
        if !isnothing(halfband)
            if (j + halfband + 1) < n
                i_range = (j+1):(j+hb)
            end
        end
        for i ∈ i_range

            for k ∈ 1:j-1
                if k == 1
                    @inbounds sumᵢ = 𝐋[i,k]*𝐋[j,k]
                else
                    @inbounds sumᵢ = sumᵢ + 𝐋[i,k]*𝐋[j,k]
                end
            end
            @inbounds 𝐋[i,j] = (𝐀[i,j] - sumᵢ) / 𝐋[j,j]
        end
    end

    #get 𝐲
    for i ∈ 1:n
        sumᵧ = 0
        for j ∈ 1:i-1
            @inbounds sumᵧ = sumᵧ + 𝐋[i,j]*𝐲[j]
        end
        @inbounds 𝐲[i] = (𝐛[i] - sumᵧ) / 𝐋[i,i]
    end
    #calculate 𝐱
    for i ∈ n:-1:1
        sumₓ = 0
        for j ∈ (i+1):n
            @inbounds sumₓ = sumₓ + 𝐋[j,i]*𝐱[j]
        end
        @inbounds 𝐱[i] = (𝐲[i] - sumₓ) / 𝐋[i,i]
    end

    return 𝐱
end

@testset "real, symmetric, and positive-definite" begin
    # test a simple 2×2 case
    test_𝐀 = [1 0; 0 1]
    test_𝐱 = [1; 1]
    test_𝐛 = test_𝐀 * test_𝐱
    @test choleski(test_𝐀, test_𝐛) ≈ test_𝐱 # approx is to account for floating point errors

    # test 10 random 2×2 cases
    for i ∈ 1:10
        test_𝐱 = rand(2)
        test_𝐛 = test_𝐀 * test_𝐱
        @test choleski(test_𝐀, test_𝐛) ≈ test_𝐱
    end

    # test a simple 3×3 case
    test_𝐀 = [2 -1 0; -1 2 -1; 0 -1 2]
    test_𝐱 = [1; 1; 1]
    test_𝐛 = [1; 0; 1]
    @test choleski(test_𝐀, test_𝐛) ≈ test_𝐱

    # test 10 random 3×3 cases
    for i ∈ 1:10
        test_𝐱 = rand(3)
        test_𝐛 = test_𝐀 * test_𝐱
        @test choleski(test_𝐀, test_𝐛) ≈ test_𝐱
    end
end