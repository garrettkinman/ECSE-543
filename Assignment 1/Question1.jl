include("Assignment 1/choleski.jl")
include("Assignment 1/circuit.jl")

# solve circuit 1 from /circuits/circuit1.toml
circuit1 = get_circuit(1)
solve_circuit(circuit1...)

# solve circuit 2 from /circuits/circuit2.toml
circuit2 = get_circuit(2)
solve_circuit(circuit2...)

# solve circuit 3 from /circuits/circuit3.toml
circuit3 = get_circuit(3)
solve_circuit(circuit3...)

# solve circuit 4 from /circuits/circuit4.toml
circuit4 = get_circuit(4)
solve_circuit(circuit4...)

# solve circuit 5 from /circuits/circuit5.toml
circuit5 = get_circuit(5)
solve_circuit(circuit5...)

# solve N×2N finite-difference mesh for N = 2
R, sJ, sR, sE = 1000, 0, 1000, 100
𝐑_equivalent = []
for N ∈ 1:32
    circuit_mesh = generate_circuit(N, R, sJ, sR, sE)
    print("N = $N\t")
    @time 𝐯_mesh =  solve_circuit(circuit_mesh...)
    𝐯₁ = 𝐯_mesh[2*N + 1]
    𝐯₂ = 0
    I = ((𝐯₂ + sE) - 𝐯₁) / sR
    R_mesh = (𝐯₁ - 𝐯₂) / I

    push!(𝐑_equivalent, R_mesh)
end