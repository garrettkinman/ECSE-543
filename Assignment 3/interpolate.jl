using Symbolics

function lagrange(𝐗, 𝐘)
    n = length(𝐗)
    @variables x # x is symbolic

    𝐋 = []

    for j ∈ 1:n
        F₁ = 1
        F₂ = 1
        for r ∈ 1:n
            # multiply in all the terms into the numerator and denominator
            if r != j
                F₁ *= (x - 𝐗[r])
                F₂ *= (𝐗[j] - 𝐗[r])
            end
        end
        
        # set the j-th Lagrange polynomial term := F (calculated above)
        push!(𝐋, F₁ / F₂)
    end
    # scale each coefficient by the corresponding y value,
    # then sum up the polynomials for the total polynomial expression
    y = sum(𝐘 .* 𝐋)

    return build_function(simplify(y, expand=true), x)
end