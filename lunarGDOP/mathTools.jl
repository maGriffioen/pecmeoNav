function diag(matrix::Array)
    sz = size(matrix)
    if sz[1] == sz[2]
        return [matrix[i, i] for i in 1:sz[1]]
    else
        error("Matrix is not square: Cannot obtain diagonal")
    end
end
norm(vector) = sqrt(sum( vector .^2 ))
