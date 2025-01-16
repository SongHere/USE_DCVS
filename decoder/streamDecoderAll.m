function [y_cell_all, M_list_all] = streamDecoderAll(GOP_list, data_stream, ...
    identifier_list_all, M_init_k, U_k, M_init_nk, U_nk, N, alpha)

[block_num, frame_num] = size(identifier_list_all);
y_cell_all = cell(block_num, frame_num);
M_list_all = zeros(block_num, frame_num);

idx_frame = 1;
idx_stream_decoding = 1;

is_skipping_list = identifier_list_all(:, idx_frame);
num_block_non_skipped = sum(is_skipping_list == 0);
y_init_stream = data_stream(idx_stream_decoding: ...
    idx_stream_decoding+num_block_non_skipped*M_init_k-1);
M_list_all(:, idx_frame) = streamDecoder(y_init_stream, ...
    identifier_list_all(:, idx_frame), M_init_k, U_k, N, alpha);
y_cell_all(:, idx_frame) = stream2Cell(data_stream(idx_stream_decoding: ...
    idx_stream_decoding+sum(M_list_all(:, idx_frame))-1), M_list_all(:, idx_frame), M_init_k);

idx_stream_decoding = idx_stream_decoding + sum(M_list_all(:, idx_frame));

idx_frame = 0;
for idx_GOP = 1:length(GOP_list)
    GOP_len = GOP_list(idx_GOP);
    decode_map = idxMidGet(1, GOP_len+1);
    decode_map = decode_map + idx_frame;

    idx_frame = decode_map(1, 3);

    num_block_non_skipped = sum(identifier_list_all(:, idx_frame) == 0);
    y_init_stream = data_stream(idx_stream_decoding: ...
        idx_stream_decoding+num_block_non_skipped*M_init_k-1);
    M_list_all(:, idx_frame) = streamDecoder(y_init_stream, ...
        identifier_list_all(:, idx_frame), M_init_k, U_k, N, alpha);
    y_cell_all(:, idx_frame) = stream2Cell(data_stream(idx_stream_decoding: ...
        idx_stream_decoding+sum(M_list_all(:, idx_frame))-1), M_list_all(:, idx_frame), M_init_k);
    idx_stream_decoding = idx_stream_decoding + sum(M_list_all(:, idx_frame));


    for idx_B_frame = 1:GOP_len - 1
        idx_frame = decode_map(idx_B_frame, 1);

        num_block_non_skipped = sum(identifier_list_all(:, idx_frame) == 0);
        y_init_stream = data_stream(idx_stream_decoding: ...
            idx_stream_decoding+num_block_non_skipped*M_init_nk-1);
        M_list_all(:, idx_frame) = streamDecoder(y_init_stream, ...
            identifier_list_all(:, idx_frame), M_init_nk, U_nk, N, alpha);
        y_cell_all(:, idx_frame) = stream2Cell(data_stream(idx_stream_decoding: ...
            idx_stream_decoding+sum(M_list_all(:, idx_frame))-1), M_list_all(:, idx_frame), M_init_nk);
        idx_stream_decoding = idx_stream_decoding + sum(M_list_all(:, idx_frame));

    end
end

end

function M_list = streamDecoder ...
    (y_init_stream, is_skipping_list, M_init, U, N, alpha)

block_num = length(is_skipping_list);
M_list = zeros(block_num, 1);

num_block_non_skipped = sum(is_skipping_list == 0);
y_init_stream = reshape(y_init_stream, [M_init, num_block_non_skipped]);
M_list_non_skipped = USE(y_init_stream, U, N, alpha);

M_list(~is_skipping_list) = M_list_non_skipped;
end
