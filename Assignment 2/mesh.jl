mutable struct PotentialMesh
    mesh::AbstractMatrix{AbstractFloat}
    h::AbstractFloat
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
    function PotentialMesh(;h::AbstractFloat=0.02)
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
        return new(mesh, h, (inner_rows, inner_cols), params...)
    end
end