% Multiple Hypothesis and Residual Reconstrution for non-skipped blocks

% references
% C. Chen, E. W. Tramel, and J. E. Fowler,“Compressed-sensing recovery of images and video using multihypothesis predictions,”in 2011 Conf. Rec. 45th Asilomar Conf. Signals, Systems and Computers (ASILOMAR), Nov. 2011, pp. 1193–1198. doi: 10.1109/ACSSC.2011.6190204.
% S. Mun and J. E. Fowler,“Residual Reconstruction for Block-Based Compressed Sensing of Video,” in 2011 Data Compression Conference, Mar. 2011, pp. 183–192. doi: 10.1109/DCC.2011.25.

function recons_block = MH_Residual_recons_IBI(y_frame, Phi_full, Psi, lambda, win_radius, ref_frame_left, ref_frame_right)

N = size(Phi_full, 2);
block_size = sqrt(N);
block_num = size(y_frame, 1);

A_full = Phi_full * Psi^(-1);

[img_rows, ~] = size(ref_frame_left);
block_row = img_rows / block_size;

win_size = 2 * win_radius + 1;

ref_frame_left = padarray(ref_frame_left, win_radius*[1, 1], 'replicate');
ref_frame_right = padarray(ref_frame_right, win_radius*[1, 1], 'replicate');

predict_block = zeros(N, block_num);
recons_block = zeros(N, block_num);

for idx_block = 1:block_num
    y = y_frame{idx_block};
    if (isempty(y))
        continue;
    end

    % MH
    % refrence window
    block_ii = rem(idx_block-1, block_row);
    block_jj = fix((idx_block - 1)/block_row);
    block_ii = block_ii * block_size + 1;
    block_jj = block_jj * block_size + 1;

    left_ref_win = ref_frame_left(block_ii:block_ii+win_size+block_size-2, ...
        block_jj:block_jj+win_size+block_size-2);
    right_ref_win = ref_frame_right(block_ii:block_ii+win_size+block_size-2, ...
        block_jj:block_jj+win_size+block_size-2);

    H_left = im2col(left_ref_win, [block_size, block_size], 'sliding');
    H_right = im2col(right_ref_win, [block_size, block_size], 'sliding');

    % Hypothesis set
    H = [H_left, H_right];

    % MH solving
    M = length(y);
    Phi = Phi_full(1:M, :);
    w = MH_ABCS(y, Phi, H, lambda);

    predict_block(:, idx_block) = H * w; % side information
    y_prediction = Phi * predict_block(:, idx_block);
    % residual reconstrution
    y_residual = y - y_prediction;

    A = A_full(1:M, :);
    recons_residual = frameReconsSPGL1(mat2cell(y_residual, M, 1), M, block_size, A);
    recons_residual = Psi^(-1) * recons_residual;
    recons_block(:, idx_block) = recons_residual + predict_block(:, idx_block);
end

end
