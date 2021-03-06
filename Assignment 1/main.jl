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

# solve NÃ2N finite-difference mesh for N = 2
R, sJ, sR, sE = 1000, 0, 1000, 100
ð_equivalent = []
for N â 1:32
    circuit_mesh = generate_circuit(N, R, sJ, sR, sE)
    halfband = 2*N + 2
    print("N = $N\tb = $halfband\t")
    @time ð¯_mesh =  solve_circuit(circuit_mesh..., halfband)
    ð¯â = ð¯_mesh[2*N + 1]
    ð¯â = 0
    I = ((ð¯â + sE) - ð¯â) / sR
    R_mesh = (ð¯â - ð¯â) / I

    push!(ð_equivalent, R_mesh)
end

## question 3

include("Assignment 1/finite-difference.jl")

# use SOR on a mesh with h=0.02 and residual limit of 1e-5
result_potentials = []
for Ï â 1.0:0.1:1.9
    mesh = PotentialMesh()
    SOR!(mesh, Ï)
    push!(result_potentials, mesh.mesh[4,3])
end
result_potentials