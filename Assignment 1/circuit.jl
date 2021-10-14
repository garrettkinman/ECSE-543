using TOML

function get_circuit(id::Int)
    circuit = TOML.tryparsefile("Assignment 1/circuits/circuit$(id).toml")
    return ( 𝐀 = circuit["A"], 𝐉 = circuit["J"], 𝐑 = circuit["R"], 𝐄 = circuit["E"])
end