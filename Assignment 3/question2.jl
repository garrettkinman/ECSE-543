using Zygote

𝐁 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
𝐇 = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.9, 8687.4, 13924.3, 22650.2]

## function declarations

function piecewise_linear(x::Real, 𝐗::AbstractVector{<:Real}, 𝐘::AbstractVector{<:Real})::Real
    n = length(𝐗)

    for i ∈ 2:n # start at second element
        if x ≤ 𝐗[i]
            slope = (𝐘[i] - 𝐘[i - 1]) / (𝐗[i] - 𝐗[i - 1]) # slop of interpolated segment
            y = slope * (x - 𝐗[i]) + 𝐘[i] # slope-intercept form
            
            return y
        elseif x > 𝐗[end] # extend out last segment
            slope = (𝐘[end] - 𝐘[end - 1]) / (𝐗[end] - 𝐗[end - 1]) # slop of last interpolated segment
            y = slope * (x - 𝐗[end]) + 𝐘[end] # slope-intercept form

            return y
        end
    end
end

f(ϕ::Real)::Real = (39.788735e6 * ϕ) + (0.3 * piecewise_linear(ϕ / 1e-4, 𝐁, 𝐇)) - 8_000

function newton_raphson(error::Real)
    # initialize estimate of flux, ϕ, as 0
    ϕₙ = 0
    i = 0

    while true
        println("Iteration: $i, Flux: $ϕₙ")

        # create new flux estimate
        ϕₙ₊₁ = ϕₙ - (f(ϕₙ) / f'(ϕₙ))
        
        # if error small enough, we're done
        # else, update the flux estimate
        if (abs(ϕₙ - ϕₙ₊₁) < error)
            break
        else
            i += 1
            ϕₙ = ϕₙ₊₁
        end
    end
end

## run newton-raphson method

newton_raphson(1e-6)