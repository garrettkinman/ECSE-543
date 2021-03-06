using Symbolics

function lagrange(š, š)
    n = length(š)
    @variables x # x is symbolic

    š = []

    for j ā 1:n
        Fā = 1
        Fā = 1
        for r ā 1:n
            # multiply in all the terms into the numerator and denominator
            if r != j
                Fā *= (x - š[r])
                Fā *= (š[j] - š[r])
            end
        end
        
        # set the j-th Lagrange polynomial term := F (calculated above)
        push!(š, Fā / Fā)
    end
    # scale each coefficient by the corresponding y value,
    # then sum up the polynomials for the total polynomial expression
    y = sum(š .* š)

    return build_function(simplify(y, expand=true), x)
end

function hermite(š, š, šā²)
    n = length(š)
    @variables x # x is symbolic

    š = []
    š = []

    for j ā 1:n
        Fā = 1
        Fā = 1
        for r ā 1:n
            # multiply in all the terms into the numerator and denominator
            if r != j
                Fā *= (x - š[r])
                Fā *= (š[j] - š[r])
            end
        end
        # calculate j-th Lagrange polynomial term, and its derivative
        L = Fā / Fā
        Lā² = Symbolics.derivative(L, x)

        # create a callable version of Lā²
        Lā±¼ā² = substitute(Lā², Dict(x => š[j]))

        push!(š, (1 - 2 * Lā±¼ā² * (x - š[j])) * (L^2))
        push!(š, (x - š[j]) * (L^2))

    end
    
    # scale each coefficient by the corresponding y or yā² value and add
    # then sum up the polynomials for the total polynomial expression
    y = sum((š .* š) + (šā² .* š))

    return build_function(simplify(y, expand=true), x)
end