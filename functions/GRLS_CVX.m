function R = GRLS_CVX(f)
%GRLS For order-1 Filter
[N, K] = size(f);

    cvx_begin
        variable R(N, N)
        minimize norm(f(:, 2:K) - R*f(:, 1:K - 1))
    cvx_end

end  

