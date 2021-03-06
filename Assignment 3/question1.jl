using Plots

include("interpolate.jl")

ð = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
ð = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]

## question 1a
@time lagrange_polynomial = lagrange(ð[1:6], ð[1:6]) |> eval

plot(0.0:0.001:1.0, lagrange_polynomial.(0.0:0.001:1.0), label="Lagrange Polynomial")
xlabel!("B (T)")
ylabel!("H (A/m)")
title!("Lagrange Polynomial H vs B")
scatter!(ð[1:6], ð[1:6], label="True Values")
savefig("Assignment 3/question1a.png")

## question 1b
@time lagrange_polynomial2 = lagrange([ð[1]; ð[9:10]; ð[13:15]], [ð[1]; ð[9:10]; ð[13:15]]) |> eval

plot(0.0:0.001:1.9, lagrange_polynomial2.(0.0:0.001:1.9), label="Lagrange Polynomial")
xlabel!("B (T)")
ylabel!("H (A/m)")
title!("Lagrange Polynomial H vs B")
scatter!(ð, ð, label="True Values")
savefig("Assignment 3/question1b.png")

## question 1c
ð_sub = [ð[1]; ð[9:10]; ð[13:15]]
ð_sub = [ð[1]; ð[9:10]; ð[13:15]]

# construct an approximate
ð_subâ² = zeros(Float64, length(ð_sub))
for j â 1:length(ð_subâ²)
    if j != length(ð_subâ²)
        ð_subâ²[j] = (ð_sub[j + 1] - ð_sub[j]) / (ð_sub[j + 1] - ð_sub[j])
    else
        ð_subâ²[j] = ð_sub[j] / ð_sub[j]
    end
end

@time hermite_polynomial = hermite(ð_sub, ð_sub, ð_subâ²) |> eval

plot(0.0:0.001:1.9, hermite_polynomial.(0.0:0.001:1.9), label="Hermite Polynomial")
xlabel!("B (T)")
ylabel!("H (A/m)")
title!("Hermite Polynomial H vs B")
scatter!(ð, ð, label="True Values")
savefig("Assignment 3/question1c.png")

##