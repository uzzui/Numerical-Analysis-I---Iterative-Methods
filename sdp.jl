include("hilbert.jl")

function s_d_p(A :: Matrix)

    is_symmetric = issymmetric(A)

    if !is_symmetric
        return false, "A matriz não é simétrica"
    end

    eigenvalues = eigvals(A)
    is_positive_definite = all(val -> val > 0, eigenvalues)

    if is_positive_definite
        return true, "A matriz é simétrica definida positiva"
    else
        return false, "A matriz não é simétrica definida positiva"
    end
end

function test_positive_definite_hilbert()
    for n in 1:100
        H, _ = hilbert(n)
        is_positive, message = s_d_p(H)

        is_positive || @info "Para n=$n: $message"
    end
end

