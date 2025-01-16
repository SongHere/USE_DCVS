function is_skipping_list = isSkipping(y_residual, M_init, tau_En_skip_mean)

En_list = sum(y_residual.^2)' / M_init;
is_skipping_list = En_list <= tau_En_skip_mean;
end