function jor(A ::Matrix, b ::Vector, x₀ ::Vector, tol ::Real, max_iter ::Integer, ω::Real)
    n = size(A, 1)
    iter = 0
    r = b - A * x₀
    r₀ = norm(r)
    err = norm(r)
    x = copy(x₀)
    
    while err > tol && iter < max_iter
        iter += 1
        
        for i in 1:n
            s = 0.0
            
            for j = 1:i-1
                s += A[i, j] * x[j]
            end
            
            for j = i+1:n
                s += A[i, j] * x[j]
            end
            
            x[i] = ω * (b[i] - s) / A[i, i] + (1 - ω) * x[i]
        end
        
        r = b - A * x
        err = norm(r) / r₀
    end
    
    return x, iter
end
