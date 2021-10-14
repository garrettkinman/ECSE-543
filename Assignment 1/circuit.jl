using TOML

"""
Retrieves a circuit of a given ID and returns its reduced incidence matrix, 𝐀; current source matrix, 𝐉;
resistance matrix, 𝐑; and voltage source matrix, 𝐄. These are returned as a named tuple.
"""
function get_circuit(id::Int)
    circuit = TOML.tryparsefile("Assignment 1/circuits/circuit$(id).toml")

    # return as named tuple
    # need to flatten 𝐀 into a matrix instead of nested vectors
    return (𝐀 = transpose(hcat(circuit["A"]...)), 𝐉 = circuit["J"], 𝐑 = circuit["R"], 𝐄 = circuit["E"])
end

"""
Solves a circuit given the reduced incidence matrix, 𝐀; current source matrix, 𝐉;
resistance matrix, 𝐑; and voltage source matrix, 𝐄. Solves using Choleski decomposition.

Tip: use the spread operator, `...`, to use the returned output of `get_circuit(id)`
without having to access the members of the named tuple.
"""
function solve_circuit(𝐀, 𝐉, 𝐑, 𝐄)
    # (𝐀𝐘𝐀ᵀ)𝐯 = 𝐀(𝐉 - 𝐘𝐄)

    # fill in 𝐘 with zeros everywhere but the main diagonal
    𝐘 = zeros(length(𝐑), length(𝐑))
    for i ∈ 1:length(𝐑)
        𝐘[i,i] = 1.0 / 𝐑[i]
    end

    # perform choleski decomposition with constructed 𝐀 and 𝐛
    𝐯 = choleski(𝐀 * 𝐘 * transpose(𝐀), 𝐀 * (𝐉 - (𝐘 * 𝐄)))
    return 𝐯
end

"""
Generates an N×2N finite-difference mesh where each branch is a resistor of given resistance.
It also has attached a branch connecting the top-right node to the bottom-left node
with a voltage source and resistor. This is to allow measurement of equivalent resistance.

Return value is in the same form as for `get_circuit`, allowing easy piping into `solve_circuit`
with the spread operator, `...`.
"""
function generate_circuit(N, R, sJ, sR, sE)
    # calculate number of branches
    n_branch = (N + 1)*2*N + N*(2*N + 1)

    # construct 𝐉, 𝐑, 𝐄 for the main N×2N mesh
    𝐉 = zeros(n_branch)
    𝐑 = R * ones(n_branch)
    𝐄 = zeros(n_branch)

    # add voltage-source branch for testing equivalent resistance across the mesh
    push!(𝐉, sJ)
    push!(𝐑, sR)
    push!(𝐄, sE)

    # calculate nodes
    n_node = (N + 1) * (2*N + 1)

    # zero-initialize 𝐀
    𝐀 = zeros(n_node, n_branch + 1)
    #print len(𝐀), len(𝐀[0])

    # set all the correct values in 𝐀
    for i_node ∈ 1:n_node
        
        # level = which row within the mesh, 1-indexed
        level = ((i_node - 1) ÷ (2*N + 1)) + 1

        # offset = which column within the mesh, 1-indexed
        offset = ((i_node - 1) % (2*N + 1)) + 1

        branch_per_level = 2*N + (2*N + 1)

        # calculate surrounding branch indices
        right = branch_per_level * (level - 1) + offset
        left = right - 1
        top = right - (2*N + 1)
        bottom = top + branch_per_level
        #print i, 'r',right, 'l', left, 't',top, 'b', bottom

        # node on top border
        if level == 1
            # top left
            if offset == 1
                𝐀[i_node, right], 𝐀[i_node, bottom] = -1, 1
            # top right
            elseif offset == (2*N + 1)
                𝐀[i_node, left], 𝐀[i_node, bottom] = 1, 1
            else
                𝐀[i_node, right], 𝐀[i_node, left], 𝐀[i_node, bottom] = -1, 1, 1
            end
        # node on bottom border
        elseif level == (N + 1)
            # bottom left
            if offset == 1
                𝐀[i_node, right], 𝐀[i_node, top] = -1, -1
            # bottom right
            elseif offset == (2*N + 1)
                𝐀[i_node, top], 𝐀[i_node, left] = -1, 1
            else
                𝐀[i_node, right], 𝐀[i_node, top], 𝐀[i_node, left] = -1, -1, 1
            end
        # node on left border
        elseif offset == 1
            𝐀[i_node, right], 𝐀[i_node, top], 𝐀[i_node, bottom] = -1, -1, 1
        # node on right border
        elseif offset == (2*N + 1)
            𝐀[i_node, top], 𝐀[i_node, left], 𝐀[i_node, bottom] = -1, 1, 1
        # nodes in middle
        else
            𝐀[i_node, right], 𝐀[i_node, top], 𝐀[i_node, left], 𝐀[i_node, bottom] = -1, -1, 1, 1
        end

        # setup the voltage source
        # right top corner
        if level == 1 && offset == (2*N + 1)
            𝐀[i_node, n_branch + 1] = -1
        elseif level == (N + 1) && offset == 1
            𝐀[i_node, n_branch + 1] = 1
        end
    end

    # set up ground node by removing that node from 𝐀
    # i.e., remove the row from 𝐀 corresponding to the bottom-left corner node
    𝐀 = 𝐀[setdiff(1:end, N*(2*N + 1) + 1), :]

    return (𝐀 = 𝐀, 𝐉 = 𝐉, 𝐑 = 𝐑, 𝐄 = 𝐄)
end