using Plots

function gauss_legendre(func::Function, interval::Tuple{Real, Real}, N::Integer; even_spacing=true, width_ratio=1)
    integral = 0

    if even_spacing
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
    else
        𝐰_relative = [width_ratio^i for i ∈ 1:N]
        𝐰 = ((interval[2] - interval[1]) / sum(𝐰_relative)) * 𝐰_relative

        𝐱 = zeros(N)
        for i ∈ 1:N
            𝐱[i] = i == 1 ? 𝐰[i] / 2 : 𝐱[i - 1] + (𝐰[i - 1] / 2) + 𝐰[i] / 2
        end

        integral = sum(𝐰 .* func.(𝐱))
    end
    
    return integral
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

ground_truth = (1*log(1) - 1)
println("Ground truth ₀∫¹ln(x)dx = $ground_truth")

println("N,Integral,Error")
for N ∈ 10:10:200
    integral = gauss_legendre(log, (0, 1), N)
    error = integral - ground_truth
    println("$N,$integral,$error")
end

## question 4c

plot(log, label="ln(x)")
plot!(cos, label="cos(x)")
xlabel!("x")
ylabel!("y")
savefig("Assignment 3/question4c.png")

ground_truth = (1*log(1) - 1)
println("Ground truth ₀∫¹ln(x)dx = $ground_truth")

println("Ratio,Integral,Error")
for r ∈ 1:0.05:2
    integral = gauss_legendre(log, (0, 1), 10, even_spacing=false, width_ratio=r)
    error = integral - ground_truth
    println("$r,$integral,$error")
end

