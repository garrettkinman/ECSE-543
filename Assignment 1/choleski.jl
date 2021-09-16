using Test

"""
Uses Choleski decomposition to solve `𝐀𝐱 = 𝐛`, where 𝐀 is real, symmetric, and positive-definite. Returns the vector 𝐱.
"""
function choleski(𝐀::AbstractMatrix{T}, 𝐛::AbstractVector{T}) where {T<:Real}
    # TODO: choleski decomposition
    # for now, use built-in solution to Ax=b
    return 𝐀 \ 𝐛
end

@testset "real, symmetric, and positive-definite" begin
    @test choleski([1 0; 0 1], [1; 1]) == [1; 1]
    @test choleski([2 -1 0; -1 2 -1; 0 -1 2], [1; 0; 1]) ≈ [1; 1; 1] # approx is to account for finit-precision floats
end