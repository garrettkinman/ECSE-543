using Test

"""
Uses Choleski decomposition to solve `๐๐ฑ = ๐`, where ๐ is real, symmetric, and positive-definite. Returns the vector ๐ฑ.
"""
function choleski(๐::AbstractMatrix{<:Real}, ๐::AbstractVector{<:Real}, halfband=nothing)
    n, m = size(๐)

    # simple error checking
    if n != m
        throw(DimensionMismatch("Matrix ๐ must be square"))
    elseif length(๐) != n
        throw(DimensionMismatch("Matrix ๐ must be nรn, and vector ๐ must be nร1"))
    end
    
    # initialize
    ๐ = zeros(n,n)
    ๐ฒ = zeros(n,1)
    ๐ฑ = zeros(n,1)

    #decompose ๐
    for j โ 1:n
        sumโฑผ = 0
        for q = 1:j-1
            @inbounds sumโฑผ = sumโฑผ + ๐[j,q]^2
        end
        @inbounds ๐[j,j] = sqrt(๐[j,j] - sumโฑผ)

        sumแตข = 0

        i_range = (j+1):n
        if !isnothing(halfband)
            if (j + halfband + 1) < n
                i_range = (j+1):(j+halfband)
            end
        end
        for i โ i_range

            for k โ 1:j-1
                if k == 1
                    @inbounds sumแตข = ๐[i,k]*๐[j,k]
                else
                    @inbounds sumแตข = sumแตข + ๐[i,k]*๐[j,k]
                end
            end
            @inbounds ๐[i,j] = (๐[i,j] - sumแตข) / ๐[j,j]
        end
    end

    #get ๐ฒ
    for i โ 1:n
        sumแตง = 0
        for j โ 1:i-1
            @inbounds sumแตง = sumแตง + ๐[i,j]*๐ฒ[j]
        end
        @inbounds ๐ฒ[i] = (๐[i] - sumแตง) / ๐[i,i]
    end
    #calculate ๐ฑ
    for i โ n:-1:1
        sumโ = 0
        for j โ (i+1):n
            @inbounds sumโ = sumโ + ๐[j,i]*๐ฑ[j]
        end
        @inbounds ๐ฑ[i] = (๐ฒ[i] - sumโ) / ๐[i,i]
    end

    return ๐ฑ
end

@testset "real, symmetric, and positive-definite" begin
    # test a simple 2ร2 case
    test_๐ = [1 0; 0 1]
    test_๐ฑ = [1; 1]
    test_๐ = test_๐ * test_๐ฑ
    @test choleski(test_๐, test_๐) โ test_๐ฑ # approx is to account for floating point errors

    # test 10 random 2ร2 cases
    for i โ 1:10
        test_๐ฑ = rand(2)
        test_๐ = test_๐ * test_๐ฑ
        @test choleski(test_๐, test_๐) โ test_๐ฑ
    end

    # test a simple 3ร3 case
    test_๐ = [2 -1 0; -1 2 -1; 0 -1 2]
    test_๐ฑ = [1; 1; 1]
    test_๐ = [1; 0; 1]
    @test choleski(test_๐, test_๐) โ test_๐ฑ

    # test 10 random 3ร3 cases
    for i โ 1:10
        test_๐ฑ = rand(3)
        test_๐ = test_๐ * test_๐ฑ
        @test choleski(test_๐, test_๐) โ test_๐ฑ
    end
end