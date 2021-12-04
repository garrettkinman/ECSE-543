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

function hermite(𝐗, 𝐘, 𝐘′)
    n = length(𝐗)
    @variables x # x is symbolic

    𝐔 = []
    𝐕 = []

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
        # calculate j-th Lagrange polynomial term, and its derivative
        L = F₁ / F₂
        L′ = Symbolics.derivative(L, x)

        # create a callable version of L′
        Lⱼ′ = substitute(L′, Dict(x => 𝐗[j]))

        push!(𝐔, (1 - 2 * Lⱼ′ * (x - 𝐗[j])) * (L^2))
        push!(𝐕, (x - 𝐗[j]) * (L^2))

    end
    
    # scale each coefficient by the corresponding y or y′ value and add
    # then sum up the polynomials for the total polynomial expression
    y = sum((𝐘 .* 𝐔) + (𝐘′ .* 𝐕))

    return build_function(simplify(y, expand=true), x)
end