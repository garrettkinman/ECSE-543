using Test
using LinearAlgebra
include("Assignment 2/mesh.jl")
include("Assignment 2/choleski.jl")

## matrix generation

"""
Generates the ๐ matrix and ๐ vector for the free nodes of a given mesh.
"""
function generate_matrix(potentials::PotentialMesh, num_nodes::Integer)
    n, m = size(potentials.mesh)

    # initialize ๐ as diagonal matrix with the -4ฯแตขโฑผ set as all the diagonals
    # need to fill in all the other 1s and 2s below
    ๐ = one(zeros(num_nodes, num_nodes)) * -4

    # zero-initialize ๐; need to fill in some of the values below
    ๐ = zeros(num_nodes)

    # k is used for indexing into our constructed ๐ and ๐
    k = 1

    # manually fill in rest of ๐ and ๐
    for i โ 1:(n-1) # right boundary are all non-free nodes
        for j โ 1:(m-1) # top boundary are all non-free nodes
            if j > 2 && potentials.mesh[i, j] == 0 && potentials.mesh[i, j - 1] == potentials.inner_potential # nodes 2, 8
                if i == 1 # node 2
                    ๐[k, k + 1] = 1
                    ๐[k, k + 2] = 2
                    ๐[k] = -potentials.inner_potential
                elseif i == 2 # node 8
                    ๐[k, k + 1] = 1
                    ๐[k, k - 2] = 1
                    ๐[k, k + 5] = 1
                    ๐[k] = -potentials.inner_potential
                end
                k += 1
            elseif j == m - 1 # nodes 3, 9, 15, 21, 27
                if i == 1 # node 3
                    ๐[k, k - 1] = 1
                    ๐[k, k + 2] = 2
                    ๐[k] = -potentials.outer_potential
                elseif i == 2 # node 9
                    ๐[k, k - 1] = 1
                    ๐[k, k + 5] = 1
                    ๐[k, k - 2] = 1
                    ๐[k] = -potentials.outer_potential
                elseif i == n - 1 # node 27
                    ๐[k, k - 1] = 1
                    ๐[k, k - 5] = 1
                    ๐[k] = -potentials.outer_potential * 2
                else # nodes 15, 21
                    ๐[k, k - 1] = 1
                    ๐[k, k + 5] = 1
                    ๐[k, k - 5] = 1
                    ๐[k] = -potentials.outer_potential
                end
                k += 1
            elseif j == 1 && i > 2 # nodes 11, 17, 23
                if potentials.mesh[i - 1, j] == potentials.inner_potential # node 11
                    ๐[k, k + 1] = 2
                    ๐[k, k + 5] = 1
                    ๐[k] = -potentials.inner_potential
                elseif i == n - 1 # node 23
                    ๐[k, k + 1] = 2
                    ๐[k, k - 5] = 1
                    ๐[k] = -potentials.outer_potential
                else # node 17
                    ๐[k, k + 1] = 2
                    ๐[k, k + 5] = 1
                    ๐[k, k - 5] = 1
                    ๐[k] = 0
                end
                k += 1
            elseif i == 3 && potentials.mesh[i - 1, j] == potentials.inner_potential # nodes 12, 13
                ๐[k, k - 1] = 1
                ๐[k, k + 1] = 1
                ๐[k, k + 5] = 1
                ๐[k] = -potentials.inner_potential
                k += 1
            elseif i == n - 1 # nodes 24, 25, 26
                ๐[k, k - 1] = 1
                ๐[k, k + 1] = 1
                ๐[k, k - 5] = 1
                ๐[k] = -potentials.outer_potential
                k += 1
            elseif i > 2 && j > 1 # nodes 14, 18, 19, 20
                ๐[k, k - 1] = 1
                ๐[k, k + 1] = 1
                ๐[k, k - 5] = 1
                ๐[k, k + 5] = 1
                ๐[k] = 0
                k += 1
            end
        end
    end
    return ๐, ๐
end

## conjugate gradient

"""
Uses Conjugate Gradient Method to solve `๐๐ฑ = ๐`, where ๐ is real and symmetric, but not necessarily positive-definite. Returns the vector ๐ฑ.
"""
function conjugate_gradient(๐::AbstractMatrix{<:Real}, ๐::AbstractVector{<:Real})
    n, m = size(๐)

    # simple error checking
    if n != m
        throw(DimensionMismatch("Matrix ๐ must be square"))
    elseif length(๐) != n
        throw(DimensionMismatch("Matrix ๐ must be nรn, and vector ๐ must be nร1"))
    end
    
    # initialize
    ๐ฑ = zeros(n,1)
    ๐ซ = ๐ - (๐ * ๐ฑ)
    ๐ฉ = copy(๐ซ)

    inf_norm_ini = 0
    two_norm_ini = 0

    for i โ 1:n
        if abs(๐ซ[i, 1]) > inf_norm_ini
            inf_norm_ini = abs(๐ซ[i, 1])
        end
        two_norm_ini += ๐ซ[i, 1]^2
    end
    two_norm_ini = โ(two_norm_ini)

    println("iteration,infinity_norm,two_norm")
    println("0,$inf_norm_ini,$two_norm_ini")

    # perform conjugate gradient method
    for i โ 1:n
        ฮฑ = (transpose(๐ฉ) * ๐ซ)[1, 1] / (transpose(๐ฉ) * ๐ * ๐ฉ)[1, 1]
        ๐ฑ = ๐ฑ + (ฮฑ * ๐ฉ)
        ๐ซ = ๐ - (๐ * ๐ฑ)
        ฮฒ = -(transpose(๐ฉ) * ๐ * ๐ซ)[1, 1] / (transpose(๐ฉ) * ๐ * ๐ฉ)[1, 1]
        ๐ฉ = ๐ซ + (ฮฒ * ๐ฉ)

        # finding the norms
        inf_norm = 0
        two_norm = 0

        for j โ 1:n
            if abs(๐ซ[j, 1]) > inf_norm
                inf_norm = abs(๐ซ[j, 1])
            end
            two_norm += ๐ซ[j, 1]^2
        end

        two_norm = โ(two_norm)
        println("$i,$inf_norm,$two_norm")
    end

    return ๐ฑ
end

## test

potentials = PotentialMesh(h=0.02)
mesh_๐, mesh_๐ = generate_matrix(potentials, 19)

conjugate_gradient(mesh_๐, mesh_๐)
conjugate_gradient(transpose(mesh_๐) * mesh_๐, transpose(mesh_๐) * mesh_๐)
choleski(transpose(mesh_๐) * mesh_๐, transpose(mesh_๐) * mesh_๐)

## test

@testset "generated ๐ and ๐ from mesh" begin
    potentials = PotentialMesh(h=0.02)
    test_๐, test_๐ = generate_matrix(potentials, 19)

    @test conjugate_gradient(transpose(test_๐) * test_๐, transpose(test_๐) * test_๐) โ choleski(transpose(test_๐) * test_๐, transpose(test_๐) * test_๐)
end