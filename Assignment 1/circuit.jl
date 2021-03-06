using TOML

"""
Retrieves a circuit of a given ID and returns its reduced incidence matrix, ๐; current source matrix, ๐;
resistance matrix, ๐; and voltage source matrix, ๐. These are returned as a named tuple.
"""
function get_circuit(id::Int)
    circuit = TOML.tryparsefile("Assignment 1/circuits/circuit$(id).toml")

    # return as named tuple
    # need to flatten ๐ into a matrix instead of nested vectors
    return (๐ = transpose(hcat(circuit["A"]...)), ๐ = circuit["J"], ๐ = circuit["R"], ๐ = circuit["E"])
end

"""
Solves a circuit given the reduced incidence matrix, ๐; current source matrix, ๐;
resistance matrix, ๐; and voltage source matrix, ๐. Solves using Choleski decomposition.

Tip: use the spread operator, `...`, to use the returned output of `get_circuit(id)`
without having to access the members of the named tuple.
"""
function solve_circuit(๐, ๐, ๐, ๐, halfband=nothing)
    # (๐๐๐แต)๐ฏ = ๐(๐ - ๐๐)

    # fill in ๐ with zeros everywhere but the main diagonal
    ๐ = zeros(length(๐), length(๐))
    for i โ 1:length(๐)
        ๐[i,i] = 1.0 / ๐[i]
    end

    # perform choleski decomposition with constructed ๐ and ๐
    ๐ฏ = choleski(๐ * ๐ * transpose(๐), ๐ * (๐ - (๐ * ๐)), halfband)
    return ๐ฏ
end

"""
Generates an Nร2N finite-difference mesh where each branch is a resistor of given resistance.
It also has attached a branch connecting the top-right node to the bottom-left node
with a voltage source and resistor. This is to allow measurement of equivalent resistance.

Return value is in the same form as for `get_circuit`, allowing easy piping into `solve_circuit`
with the spread operator, `...`.
"""
function generate_circuit(N, R, sJ, sR, sE)
    # calculate number of branches
    n_branch = (N + 1)*2*N + N*(2*N + 1)

    # construct ๐, ๐, ๐ for the main Nร2N mesh
    ๐ = zeros(n_branch)
    ๐ = R * ones(n_branch)
    ๐ = zeros(n_branch)

    # add voltage-source branch for testing equivalent resistance across the mesh
    push!(๐, sJ)
    push!(๐, sR)
    push!(๐, sE)

    # calculate nodes
    n_node = (N + 1) * (2*N + 1)

    # zero-initialize ๐
    ๐ = zeros(n_node, n_branch + 1)
    #print len(๐), len(๐[0])

    # set all the correct values in ๐
    for i_node โ 1:n_node
        
        # level = which row within the mesh, 1-indexed
        level = ((i_node - 1) รท (2*N + 1)) + 1

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
                ๐[i_node, right], ๐[i_node, bottom] = -1, 1
            # top right
            elseif offset == (2*N + 1)
                ๐[i_node, left], ๐[i_node, bottom] = 1, 1
            else
                ๐[i_node, right], ๐[i_node, left], ๐[i_node, bottom] = -1, 1, 1
            end
        # node on bottom border
        elseif level == (N + 1)
            # bottom left
            if offset == 1
                ๐[i_node, right], ๐[i_node, top] = -1, -1
            # bottom right
            elseif offset == (2*N + 1)
                ๐[i_node, top], ๐[i_node, left] = -1, 1
            else
                ๐[i_node, right], ๐[i_node, top], ๐[i_node, left] = -1, -1, 1
            end
        # node on left border
        elseif offset == 1
            ๐[i_node, right], ๐[i_node, top], ๐[i_node, bottom] = -1, -1, 1
        # node on right border
        elseif offset == (2*N + 1)
            ๐[i_node, top], ๐[i_node, left], ๐[i_node, bottom] = -1, 1, 1
        # nodes in middle
        else
            ๐[i_node, right], ๐[i_node, top], ๐[i_node, left], ๐[i_node, bottom] = -1, -1, 1, 1
        end

        # setup the voltage source
        # right top corner
        if level == 1 && offset == (2*N + 1)
            ๐[i_node, n_branch + 1] = -1
        elseif level == (N + 1) && offset == 1
            ๐[i_node, n_branch + 1] = 1
        end
    end

    # set up ground node by removing that node from ๐
    # i.e., remove the row from ๐ corresponding to the bottom-left corner node
    ๐ = ๐[setdiff(1:end, N*(2*N + 1) + 1), :]

    return (๐ = ๐, ๐ = ๐, ๐ = ๐, ๐ = ๐)
end