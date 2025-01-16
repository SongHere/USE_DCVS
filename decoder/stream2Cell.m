function y_cell = stream2Cell(y_stream, M_list, M_init)

block_num = length(M_list);

y_cell = cell(block_num, 1);
y_init = zeros(M_init, block_num);

is_non_skipped_list = M_list > 0;
num_block_non_skipped = sum(is_non_skipped_list);

y_init(:, is_non_skipped_list) = reshape(y_stream(1:num_block_non_skipped*M_init), [M_init, num_block_non_skipped]);
idx_decoding = num_block_non_skipped * M_init + 1;

for idx_block = 1:block_num
    M = M_list(idx_block);
    if (M == 0)
        continue;
    end
    M_add = M - M_init;
    y_cell{idx_block} = [y_init(:,idx_block); y_stream(idx_decoding:idx_decoding+M_add-1)];
    idx_decoding = idx_decoding+M_add;

end

end