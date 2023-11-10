function cholesky_decomposition(A :: Matrix)
    n = size(A, 1)
    L = zeros(BigFloat, n, n)

    for i in 1:n
        for j in 1:i
            sum_val = 0.0
            for k in 1:j-1
                sum_val += L[i, k] * L[j, k]
            end

            if i == j
                L[i, j] = sqrt(A[i, i] - sum_val)
            else
                L[i, j] = (1.0 / L[j, j]) * (A[i, j] - sum_val)
            end
        end
    end

    return L
end

function solve_Cholesky(A :: Matrix ,b :: Vector)
	L = cholesky_decomposition(A)
	y = L \ b
	x = transpose(L) \ y
	return x
end
