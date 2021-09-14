using Test

"""
Uses Choleski decomposition to solve `𝐀𝐱 = 𝐛`, where 𝐀 is real, symmetric, and positive-definite. Returns the vector 𝐱.
"""
function choleski(𝐀::AbstractMatrix{Real}, 𝐛::AbstractVector{Real})::AbstractVector{Real}
    # TODO: choleski decomposition
    # for now, use built-in solution to Ax=b
    return 𝐀 \ 𝐛
end

