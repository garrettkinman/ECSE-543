const h = 0.02
const MESH_NODES = [0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 0 0 0 0
                    0 0 0 15 15 15
                    0 0 0 15 15 15]

function SOR(ω::AbstractFloat)
    # we should keep in mind that the most 'outer' node has V = 0
    for y ∈ 1:numRows
        for x ∈ 1:numColumns
            if x == 0 && y > innerHeight / h
                # we can assume this formula since when x = 0 d(potential)/dx = 0, we have mesh[x - 1][y] = mesh[x+1}[y]
                MESH_NODES[y, x] = (1 - ω) * MESH_NODES[y, x] + (ω / 4) * (2*MESH_NODES[y, x+1] + MESH_NODES[y-1, x] + MESH_NODES[y+1, x])
            elseif y == 0 && x > innerWidth / h
                # same argument as above
                MESH_NODES[y, x] = (1 - ω) * MESH_NODES[y, x] + (ω / 4) * (MESH_NODES[y, x-1] + MESH_NODES[y, x+1] + 2*MESH_NODES[y+1, x])
            elseif x > innerWidth / h || y > innerHeight / h
                MESH_NODES[y, x] = (1 - ω) * MESH_NODES[y, x] + (ω / 4) * (MESH_NODES[y, x-1] + MESH_NODES[y, x+1] + MESH_NODES[y-1, x] + MESH_NODES[y+1, x])
            end
        end
    end
    return MESH_NODES
end