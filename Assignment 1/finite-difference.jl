mutable struct PotentialMesh
    mesh::AbstractMatrix{AbstractFloat}
    h::AbstractFloat
    residual_limit::AbstractFloat
    inner_size::Tuple{Integer, Integer}
    outer_height::AbstractFloat
    outer_width::AbstractFloat
    outer_potential::AbstractFloat
    inner_height::AbstractFloat
    inner_width::AbstractFloat
    inner_potential::AbstractFloat

    """
    Constructor for PotentialMesh with optional kwargs for `h`, the mesh spacing, and `residual_limit`,
    the cutoff point for SOR.
    """
    function PotentialMesh(;h::AbstractFloat=0.02, residual_limit::AbstractFloat=1e-5)
        # constant parameters for defining the operating geometry and potentials
        params = (outer_height=0.1, outer_width=0.1, outer_potential=0.0, inner_height=0.02, inner_width=0.04, inner_potential=15.0)

        # construct a matrix where each element is a node at the intersection of neighboring meshes
        # initialize values to outer_potential
        num_rows, num_cols = Int.(((params.outer_height / h) + 1, (params.outer_width / h) + 1))
        mesh = fill(params.outer_potential, (num_rows, num_cols))

        # fill in the nodes corresponding to the inner conductor with inner_potential
        inner_rows, inner_cols = Int.(((params.inner_height / h) + 1, (params.inner_width / h) + 1))
        mesh[1:inner_rows, 1:inner_cols] .= params.inner_potential

        # construct full struct with the calculated and constant parameters
        return new(mesh, h, residual_limit, (inner_rows, inner_cols), params...)
    end
end

"""
Performs one iteration of SOR on a potential mesh.
"""
function SOR_iterate!(potentials::PotentialMesh, ω::AbstractFloat)
    num_rows, num_cols = size(potentials.mesh)
    inner_rows, inner_cols = potentials.inner_size

    for j ∈ 1:num_rows-1
        for i ∈ 1:num_cols-1
            if i == 1 && j > inner_rows # if along y-axis but not in inner conductor
                # because of symmetry, we can assume the potential at [j,i-1] = [j,i+1], hence the 2*[j,i+1]
                potentials.mesh[j, i] = (1 - ω) * potentials.mesh[j, i] + (ω / 4) * (2*potentials.mesh[j, i+1] + potentials.mesh[j-1, i] + potentials.mesh[j+1, i])
            elseif j == 1 && i > inner_cols # if along x-axis but not in inner conductor
                # same argument as above, just changed for x-axis
                potentials.mesh[j, i] = (1 - ω) * potentials.mesh[j, i] + (ω / 4) * (potentials.mesh[j, i-1] + potentials.mesh[j, i+1] + 2*potentials.mesh[j+1, i])
            elseif i > inner_cols || j > inner_rows # if somewhere else outside inner conductor
                potentials.mesh[j, i] = (1 - ω) * potentials.mesh[j, i] + (ω / 4) * (potentials.mesh[j, i-1] + potentials.mesh[j, i+1] + potentials.mesh[j-1, i] + potentials.mesh[j+1, i])
            end # don't do anything for anything in the inner conductor
        end
    end

    return potentials
end

"""
Calculates the residual of a given potential mesh state.
"""
function get_residual(potentials::PotentialMesh)
    residual_cur = 0
    residual_fin = 0

    num_rows, num_cols = size(potentials.mesh)
    inner_rows, inner_cols = potentials.inner_size

    for j ∈ 1:num_rows-1
        for i ∈ 1:num_cols-1
            if i == 1 && j > inner_rows # if along y-axis but not in inner conductor
                # because of symmetry, we can assume the potential at [j,i-1] = [j,i+1], hence the 2*[j,i+1]
                residual_cur = 2*potentials.mesh[j, i+1] + potentials.mesh[j-1, i] + potentials.mesh[j+1, i] - 4*potentials.mesh[j, i]
            elseif j == 1 && i > inner_cols # if along x-axis but not in inner conductor
                # same argument as above, just changed for x-axis
                residual_cur = potentials.mesh[j, i-1] + potentials.mesh[j, i+1] + 2*potentials.mesh[j+1, i] - 4*potentials.mesh[j, i]
            elseif i > inner_cols || j > inner_rows # if somewhere else outside inner conductor
                residual_cur = potentials.mesh[j, i-1] + potentials.mesh[j, i+1] + potentials.mesh[j-1, i] + potentials.mesh[j+1, i] - 4*potentials.mesh[j, i]
            end # don't do anything for anything in the inner conductor

            residual_cur = abs(residual_cur)
            if residual_cur > residual_fin
                # update variable with the biggest residue amongst the free point
                residual_fin = residual_cur
            end
        end
    end
    return residual_fin
end

"""
Performs the whole SOR algorithm.
"""
function SOR!(potentials::PotentialMesh, ω::AbstractFloat)
    iteration = 0
    SOR_iterate!(potentials, ω)
    while get_residual(potentials) >= potentials.residual_limit
        SOR_iterate!(potentials, ω)
        iteration += 1
    end
    println("total iteration is: $iteration")
    return potentials
end