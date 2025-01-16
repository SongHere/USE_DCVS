function y_stream_add = supplementarySampling(img_block_non_skipped, Phi_add, M_add)
% resampling
block_num = length(M_add);
y_stream_add = zeros(sum(M_add), 1); % measurements stream

idx_samples = 1;
for idx_block = 1:block_num
    M = M_add(idx_block);
    y_add = Phi_add(1:M, :) * img_block_non_skipped(:, idx_block);
    y_stream_add(idx_samples: idx_samples+M-1) = y_add;
    idx_samples = idx_samples + M;
end

end