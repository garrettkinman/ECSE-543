using TOML

"""
Retrieves a circuit of a given ID and returns its reduced incidence matrix, 𝐀; current source matrix, 𝐉;
resistance matrix, 𝐑; and voltage source matrix, 𝐄. These are returned as a named tuple.
"""
function get_circuit(id::Int)
    circuit = TOML.tryparsefile("Assignment 1/circuits/circuit$(id).toml")

    # return as named tuple
    # need to flatten 𝐀 into a matrix instead of nested vectors
    return ( 𝐀 = transpose(hcat(circuit["A"]...)), 𝐉 = circuit["J"], 𝐑 = circuit["R"], 𝐄 = circuit["E"])
end

"""
Solves a circuit given the reduced incidence matrix, 𝐀; current source matrix, 𝐉;
resistance matrix, 𝐑; and voltage source matrix, 𝐄. Solves using Choleski decomposition.

Tip: use the spread operator, `...`, to use the returned output of `get_circuit(id)`
without having to access the members of the named tuple.
"""
function solve_circuit(𝐀, 𝐉, 𝐑, 𝐄)
    # (𝐀 * 𝐘 * 𝐀ᵀ) * 𝐯 = 𝐀 * (𝐉 - 𝐘 * 𝐄)

    # fill in 𝐘 with zeros everywhere but the main diagonal
    𝐘 = zeros(length(𝐑), length(𝐑))
    for i ∈ 1:length(𝐑)
        𝐘[i,i] = 1.0 / 𝐑[i]
    end

    # perform choleski decomposition with constructed 𝐀 and 𝐛
    𝐯 = choleski(𝐀 * 𝐘 * transpose(𝐀), 𝐀 * (𝐉 - (𝐘 * 𝐄)))
    return 𝐯
end