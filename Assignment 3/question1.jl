using Plots

include("interpolate.jl")

𝐁 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
𝐇 = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]

## question 1a
@time lagrange_polynomial = lagrange(𝐁[1:6], 𝐇[1:6]) |> eval

plot(0.0:0.001:1.0, lagrange_polynomial.(0.0:0.001:1.0), label="Lagrange Polynomial")
xlabel!("B (T)")
ylabel!("H (A/m)")
title!("Lagrange Polynomial H vs B")
scatter!(𝐁[1:6], 𝐇[1:6], label="True Values")
savefig("Assignment 3/question1a.png")

## question 1b
@time lagrange_polynomial2 = lagrange([𝐁[1]; 𝐁[9:15]], [𝐇[1]; 𝐇[9:15]]) |> eval

plot(0.0:0.001:1.9, lagrange_polynomial2.(0.0:0.001:1.9), label="Lagrange Polynomial")
xlabel!("B (T)")
ylabel!("H (A/m)")
title!("Lagrange Polynomial H vs B")
scatter!([𝐁[1]; 𝐁[9:15]], [𝐇[1]; 𝐇[9:15]], label="True Values")
savefig("Assignment 3/question1b.png")
