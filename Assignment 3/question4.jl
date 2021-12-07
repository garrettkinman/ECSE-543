
function gauss_legendre(func::Function, interval::Tuple{Real, Real}, N::Integer)
    segment_width = (interval[2] - interval[1]) / N

    # each element, xᵢ, of 𝐱 will hold the midpoint of the i-th segment
    # e.g., a segment from x = 0.1 to 0.2 will be stored as 0.1
    𝐱 = zeros(N)
    for i ∈ 1:N
        𝐱[i] = (segment_width * i) - (segment_width / 2)
    end

    # each element, wᵢ, of 𝐰 will hold the width of the i-th segment
    𝐰 = ones(N) * segment_width

    integral = sum(𝐰 .* func.(𝐱))
end

## question 4a

ground_truth = sin(1) - sin(0)
println("Ground truth ₀∫¹cos(x)dx = $ground_truth")

println("N,Integral,Error")
for N ∈ 1:20
    integral = gauss_legendre(cos, (0, 1), N)
    error = integral - ground_truth
    println("$N,$integral,$error")
end

## question 4b

ground_truth = sin(1) - sin(0)
println("Ground truth ₀∫¹cos(x)dx = $ground_truth")

println("N,Integral,Error")
for N ∈ 10:10:200
    integral = gauss_legendre(cos, (0, 1), N)
    error = integral - ground_truth
    println("$N,$integral,$error")
end