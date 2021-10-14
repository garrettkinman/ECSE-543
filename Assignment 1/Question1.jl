include("Assignment 1/choleski.jl")
include("Assignment 1/circuit.jl")

# solve circuit 1 from /circuits/circuit1.toml
circuit = get_circuit(1)
solve_circuit(circuit...)