using Zygote

f₁(v₁::Real, v₂::Real) = 44.06e-5 - (0.002 * v₁) - (6e-7 * exp(40 * (v₁ - v₂)))
f₂(v₁::Real, v₂::Real) = 44.12e-5 - (0.002 * v₁) - (12e-7 * exp(40 * v₂))

function newton_raphson(error::Real)
    # initialize estimate of node voltages, 𝐯, as 0
    𝐯ₙ = [0, 0]
    i = 0

    println("iteration,v₁,v₂,error")

    while true
        # create new node voltage estimates
        # 𝐯ₙ₊₁ = 𝐯ₙ - (𝐟′⁻¹·𝐟)
        𝐟′ = [
            gradient(x -> f₁(x, 𝐯ₙ[2]), 𝐯ₙ[1])[1] gradient(x -> f₁(𝐯ₙ[1], x), 𝐯ₙ[2])[1]
            gradient(x -> f₂(x, 𝐯ₙ[2]), 𝐯ₙ[1])[1] gradient(x -> f₂(𝐯ₙ[1], x), 𝐯ₙ[2])[1]
        ]

        𝐯ₙ₊₁ = 𝐯ₙ - (inv(𝐟′) * [f₁(𝐯ₙ...), f₂(𝐯ₙ...)])

        ϵₖ = max(abs.(𝐯ₙ - 𝐯ₙ₊₁)...)

        println("$i,$(𝐯ₙ[1]),$(𝐯ₙ[2]),$ϵₖ")
        
        # if error small enough, we're done
        # else, update the node voltage estimates
        if (ϵₖ < error)
            break
        else
            i += 1
            𝐯ₙ = 𝐯ₙ₊₁
        end
    end
end

newton_raphson(1e-5)