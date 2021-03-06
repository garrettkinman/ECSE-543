using Plots

function gauss_legendre(func::Function, interval::Tuple{Real, Real}, N::Integer; even_spacing=true, width_ratio=1)
    integral = 0

    if even_spacing
        segment_width = (interval[2] - interval[1]) / N

        # each element, xแตข, of ๐ฑ will hold the midpoint of the i-th segment
        # e.g., a segment from x = 0.1 to 0.2 will be stored as 0.1
        ๐ฑ = zeros(N)
        for i โ 1:N
            ๐ฑ[i] = (segment_width * i) - (segment_width / 2)
        end
    
        # each element, wแตข, of ๐ฐ will hold the width of the i-th segment
        ๐ฐ = ones(N) * segment_width
    
        integral = sum(๐ฐ .* func.(๐ฑ))
    else
        ๐ฐ_relative = [width_ratio^i for i โ 1:N]
        ๐ฐ = ((interval[2] - interval[1]) / sum(๐ฐ_relative)) * ๐ฐ_relative

        ๐ฑ = zeros(N)
        for i โ 1:N
            ๐ฑ[i] = i == 1 ? ๐ฐ[i] / 2 : ๐ฑ[i - 1] + (๐ฐ[i - 1] / 2) + ๐ฐ[i] / 2
        end

        integral = sum(๐ฐ .* func.(๐ฑ))
    end
    
    return integral
end

## question 4a

ground_truth = sin(1) - sin(0)
println("Ground truth โโซยนcos(x)dx = $ground_truth")

println("N,Integral,Error")
for N โ 1:20
    integral = gauss_legendre(cos, (0, 1), N)
    error = integral - ground_truth
    println("$N,$integral,$error")
end

## question 4b

ground_truth = (1*log(1) - 1)
println("Ground truth โโซยนln(x)dx = $ground_truth")

println("N,Integral,Error")
for N โ 10:10:200
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
println("Ground truth โโซยนln(x)dx = $ground_truth")

println("Ratio,Integral,Error")
for r โ 1:0.05:2
    integral = gauss_legendre(log, (0, 1), 10, even_spacing=false, width_ratio=r)
    error = integral - ground_truth
    println("$r,$integral,$error")
end