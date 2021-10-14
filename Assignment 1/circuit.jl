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