function M_list = USE(y_init, U, N, alpha)
% Unified Sparsity Estimation

[M_init, ~] = size(y_init);

kappa = U * y_init;
k_list = sum(abs(kappa) >= 1)' + 1;

% measurement allocation
M_list = round(alpha*(9.62 * log2(k_list) + 0.85));
M_list(M_list > N) = N;
M_list(M_list < M_init) = M_init;

end