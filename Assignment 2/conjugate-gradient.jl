using Test
using LinearAlgebra
include("Assignment 2/mesh.jl")
include("Assignment 2/choleski.jl")

function generate_matrix(potentials::PotentialMesh, num_nodes::Integer)
    n, m = size(potentials.mesh)

    # initialize 𝐀 as diagonal matrix with the -4φᵢⱼ set as all the diagonals
    # need to fill in all the other 1s and 2s below
    𝐀 = one(zeros(num_nodes, num_nodes)) * -4

    # zero-initialize 𝐛; need to fill in some of the values below
    𝐛 = zeros(num_nodes)

    # k is used for indexing into our constructed 𝐀 and 𝐛
    k = 1

    # manually fill in rest of 𝐀 and 𝐛
    for i ∈ 1:(n-1) # right boundary are all non-free nodes
        for j ∈ 1:(m-1) # top boundary are all non-free nodes
            if j > 2 && potentials.mesh[i, j] == 0 && potentials.mesh[i, j - 1] == potentials.inner_potential # nodes 2, 8
                if i == 1 # node 2
                    𝐀[k, k + 1] = 1
                    𝐀[k, k + 2] = 2
                    𝐛[k] = -potentials.inner_potential
                elseif i == 2 # node 8
                    𝐀[k, k + 1] = 1
                    𝐀[k, k - 2] = 1
                    𝐀[k, k + 5] = 1
                    𝐛[k] = -potentials.inner_potential
                end
                k += 1
            elseif j == m - 1 # nodes 3, 9, 15, 21, 27
                if i == 1 # node 3
                    𝐀[k, k - 1] = 1
                    𝐀[k, k + 2] = 2
                    𝐛[k] = -potentials.outer_potential
                elseif i == 2 # node 9
                    𝐀[k, k - 1] = 1
                    𝐀[k, k + 5] = 1
                    𝐀[k, k - 2] = 1
                    𝐛[k] = -potentials.outer_potential
                elseif i == n - 1 # node 27
                    𝐀[k, k - 1] = 1
                    𝐀[k, k - 5] = 1
                    𝐛[k] = -potentials.outer_potential * 2
                else # nodes 15, 21
                    𝐀[k, k - 1] = 1
                    𝐀[k, k + 5] = 1
                    𝐀[k, k - 5] = 1
                    𝐛[k] = -potentials.outer_potential
                end
                k += 1
            elseif j == 1 && i > 2 # nodes 11, 17, 23
                if potentials.mesh[i - 1, j] == potentials.inner_potential # node 11
                    𝐀[k, k + 1] = 2
                    𝐀[k, k + 5] = 1
                    𝐛[k] = -potentials.inner_potential
                elseif i == n - 1 # node 23
                    𝐀[k, k + 1] = 2
                    𝐀[k, k - 5] = 1
                    𝐛[k] = -potentials.outer_potential
                else # node 17
                    𝐀[k, k + 1] = 2
                    𝐀[k, k + 5] = 1
                    𝐀[k, k - 5] = 1
                    𝐛[k] = 0
                end
                k += 1
            elseif i == 3 && potentials.mesh[i - 1, j] == potentials.inner_potential # nodes 12, 13
                𝐀[k, k - 1] = 1
                𝐀[k, k + 1] = 1
                𝐀[k, k + 5] = 1
                𝐛[k] = -potentials.inner_potential
                k += 1
            elseif i == n - 1 # nodes 24, 25, 26
                𝐀[k, k - 1] = 1
                𝐀[k, k + 1] = 1
                𝐀[k, k - 5] = 1
                𝐛[k] = -potentials.outer_potential
                k += 1
            elseif i > 2 && j > 1 # nodes 14, 18, 19, 20
                𝐀[k, k - 1] = 1
                𝐀[k, k + 1] = 1
                𝐀[k, k - 5] = 1
                𝐀[k, k + 5] = 1
                𝐛[k] = 0
                k += 1
            end
        end
    end
    return 𝐀, 𝐛
end

## debug

potentials = PotentialMesh(h=0.02)
mesh_𝐀, mesh_𝐛 = generate_matrix(potentials, 19)

## conjugate gradient

"""
Uses Conjugate Gradient Descent to solve `𝐀𝐱 = 𝐛`, where 𝐀 is real, symmetric, and positive-definite. Returns the vector 𝐱.
"""
function conjugate_gradient(𝐀::AbstractMatrix{<:Real}, 𝐛::AbstractVector{<:Real})
    n, m = size(𝐀)

    # simple error checking
    if n != m
        throw(DimensionMismatch("Matrix 𝐀 must be square"))
    elseif length(𝐛) != n
        throw(DimensionMismatch("Matrix 𝐀 must be n×n, and vector 𝐛 must be n×1"))
    end
    
    𝐱 = zeros(n,1)
    𝐫 = 𝐛 - (𝐀 * 𝐱)
    𝐩 = copy(𝐫)

    inf_norm_ini = 0
    two_norm_ini = 0

    # debug
    println(𝐫)
    println(typeof(𝐫))
    println(size(𝐫))

    for i ∈ 1:n
        if abs(𝐫[i, 1]) > inf_norm_ini
            inf_norm_ini = abs(𝐫[i, 1])
        end
        two_norm_ini += 𝐫[i, 1]^2
    end
    two_norm_ini = √(two_norm_ini)

    println("iteration,infinity_norm,two_norm")
    println("0,$inf_norm_ini,$two_norm_ini")

    for i ∈ 1:n
        α = (transpose(𝐩) * 𝐫)[1, 1] / (transpose(𝐩) * 𝐀 * 𝐩)[1, 1]
        𝐱 = 𝐱 + (α * 𝐩)
        𝐫 = 𝐛 - (𝐀 * 𝐱)
        β = -(transpose(𝐩) * 𝐀 * 𝐫)[1, 1] / (transpose(𝐩) * 𝐀 * 𝐩)[1, 1]
        𝐩 = 𝐫 + (β * 𝐩)

        # finding the norms
        inf_norm = 0
        two_norm = 0

        for j ∈ 1:n
            if abs(𝐫[j, 1]) > inf_norm
                inf_norm = abs(𝐫[j, 1])
            end
            two_norm += 𝐫[j, 1]^2
        end

        two_norm = √(two_norm)
        println("$i,$inf_norm,$two_norm")
    end

    return 𝐱
end

## debug

# test_𝐀 = [1 0; 0 1]
# test_𝐱 = [1; 1]
# test_𝐛 = test_𝐀 * test_𝐱

# @test conjugate_gradient(test_𝐀, test_𝐛) ≈ test_𝐱

conjugate_gradient(transpose(mesh_𝐀) * mesh_𝐀, transpose(mesh_𝐀) * mesh_𝐛)
choleski(transpose(mesh_𝐀) * mesh_𝐀, transpose(mesh_𝐀) * mesh_𝐛)
## test

@testset "real, symmetric, and positive-definite" begin
    # test a simple 2×2 case
    test_𝐀 = [1 0; 0 1]
    test_𝐱 = [1; 1]
    test_𝐛 = test_𝐀 * test_𝐱
    @test conjugate_gradient(test_𝐀, test_𝐛) ≈ test_𝐱 # approx is to account for floating point errors

    # test 10 random 2×2 cases
    for i ∈ 1:10
        test_𝐱 = rand(2)
        test_𝐛 = test_𝐀 * test_𝐱
        @test conjugate_gradient(test_𝐀, test_𝐛) ≈ test_𝐱
    end

    # test a simple 3×3 case
    test_𝐀 = [2 -1 0; -1 2 -1; 0 -1 2]
    test_𝐱 = [1; 1; 1]
    test_𝐛 = [1; 0; 1]
    @test conjugate_gradient(test_𝐀, test_𝐛) ≈ test_𝐱

    # test 10 random 3×3 cases
    for i ∈ 1:10
        test_𝐱 = rand(3)
        test_𝐛 = test_𝐀 * test_𝐱
        @test conjugate_gradient(test_𝐀, test_𝐛) ≈ test_𝐱
    end
end