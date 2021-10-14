include("Assignment 1/choleski.jl")
include("Assignment 1/circuit.jl")

# solve circuit 1 from /circuits/circuit1.toml

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

circuit = get_circuit(1)

solve_circuit(circuit.𝐀, circuit.𝐉, circuit.𝐑, circuit.𝐄)