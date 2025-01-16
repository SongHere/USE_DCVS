% -----------------------------------------------------------------------
% The MATLAB programs of (decoder)
% Unified Sparsity Estimation-Based Distributed Compressive Video Sensing

% Author: Zhen Song
% The results within the paper are test on MATLAB R2018a
% -----------------------------------------------------------------------

%% 
% running the encoder first

addpath('.\decoder');
% download the 'spgl1' package and install
% https://friedlander.io/spgl1/
% addpath('.\decoder\spgl1-2.1')

%% arguments configration
% params require transmitted from the encoder

% GOP_len_set = [2, 4, 8, 16];
% GOP_len_init = GOP_len_set(2);
% tau_GOP_set = [0.1, 0.3, 0.6];

% qf_k = 60;
% qf_nk = 30;
% alpha = 1.7;
% block_size = 8;
% M_init_k = 27;
% M_init_nk = 18;
% Mc_k = 32;
% Mc_nk = 24;

% above params should be modified along with the encoder

% calculable params
N = block_size^2;
block_num = height * width / N;

%% matrix construction through random seed synchronization
% the construction methods are the same as the encoder
% random measruement matrix
Phi_state = 0; % the random seed
randn('seed', Phi_state);
Phi_full = randn(N, N);
Phi_full = orth(Phi_full)';

% JPEG 8*8 luminance quantization table
q_table = JPEG88QTable();
[~, q_table, idx_col2zig] = zigzagScanning(q_table);
% q_table: JPEG 8*8 luminance quantization table
% idx_col2zig: index relationship between raster scanning (value) and zig-zag scanning (position)

Q_k = qTableAdjustment(q_table, qf_k); % threshold adjustment
Q_nk = qTableAdjustment(q_table, qf_nk);

% DCT basis for saprse representation
load Psi_DCT88;

% initial linear reconstruction
Phi_k = Phi_full(1:M_init_k, :);
Phi_nk = Phi_full(1:M_init_nk, :);

% pixel domian correlation matrix
num_Rxx = 0.95;
Rxx = ImageCorr(block_size, block_size, num_Rxx);
% linear recons matrix construction
P_k = Rxx * Phi_k' * (Phi_k * Rxx * Phi_k')^(-1);
P_nk = Rxx * Phi_nk' * (Phi_nk * Rxx * Phi_nk')^(-1);

% Selected USE (in zig-zag scanning order)
I_N = eye(N); % identity matrix

idx_zig_k = idx_col2zig(2:Mc_k); % selected index
Lambda_k = I_N(idx_zig_k, :); % permutation matrix
Gamma_k = I_N .* Q_k; % threshold matrix
U_k = Lambda_k * (Gamma_k \ (Psi * P_k)); % unified sparsity estimation matrix

idx_zig_nk = idx_col2zig(2:Mc_nk);
Lambda_nk = I_N(idx_zig_nk, :);
Gamma_nk = I_N .* Q_nk;
U_nk = Lambda_nk * (Gamma_nk \ (Psi * P_nk));

% above matrices are used for measurement stream decoding

A_full = Phi_full * Psi^(-1); % for sparse decoding

%% measurement stream decoding
% implements the USE independently
GOP_list = GOPListDecoder(GOP_len_init, tau_GOP_set, GOP_len_set, identifier_list_all);
% GOP size determines the type of each frame
[y_cell_all, M_list_all] = streamDecoderAll(GOP_list, data_stream, ...
    identifier_list_all, M_init_k, U_k, M_init_nk, U_nk, N, alpha);

%% frame reconstruction
% reconstruction result variable
recons_block_seq = zeros(N, block_num, frame_num);
recons_seq = zeros(height, width, frame_num); % reconstruction seq result
psnr_list = zeros(frame_num, 1); % psnr of each frame
ssim_list = psnr_list;

% MH
lambda = 0.10;
MH_win_radius = 4;

%% decoding
% first frame
idx_frame = 1;

% reconstrution
recons_block_non_skipped = frameReconsSPGL1(y_cell_all(:, idx_frame), M_list_all(:, idx_frame), block_size, A_full);
recons_block_non_skipped = Psi^(-1) * recons_block_non_skipped;

% fine-tuning
recons_block_non_skipped(recons_block_non_skipped > 255) = 255;
recons_block_non_skipped(recons_block_non_skipped < 0) = 0;

recons_block_seq(:, :, idx_frame) = recons_block_non_skipped;
recons_seq(:, :, idx_frame) = col2im(recons_block_non_skipped, [block_size, block_size], [height, width], 'distinct');

psnr_list(idx_frame) = Psnr(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
ssim_list(idx_frame) = ssim(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
fprintf('GOP: 1, Size = %d\n I Frame: %d, PSNR = %.2f dB, SSIM = %.4f \n', ...
    GOP_len_init, idx_frame, psnr_list(idx_frame), ssim_list(idx_frame));

idx_frame = 0;
for idx_GOP = 1:length(GOP_list)
    GOP_len = GOP_list(idx_GOP);
    decoding_map = idxMidGet(1, GOP_len+1); % = encoding_map
    decoding_map = decoding_map + idx_frame;

    % I/P frame
    idx_frame = decoding_map(1, 3);

    % MH and Residual reconstruction
    recons_block_non_skipped = MH_Residual_recons_IBI(y_cell_all(:, idx_frame), Phi_full, Psi, ...
        lambda, MH_win_radius, recons_seq(:, :, idx_frame-GOP_len), ...
        zeros(height, width));

    is_skipping_list = identifier_list_all(:, idx_frame) > 0;
    recons_block = recons_block_non_skipped;
    recons_block(:, is_skipping_list) = recons_block_seq(:, is_skipping_list, idx_frame-GOP_len);

    recons_block(recons_block > 255) = 255;
    recons_block(recons_block < 0) = 0;
    recons_block_seq(:, :, idx_frame) = recons_block;
    recons_seq(:, :, idx_frame) = col2im(recons_block, [block_size, block_size], [height, width], 'distinct');

    psnr_list(idx_frame) = Psnr(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
    ssim_list(idx_frame) = ssim(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
    fprintf('GOP: %d, Size = %d\n I/P Frame: %d, PSNR = %.2f dB, SSIM = %.4f \n', ...
        idx_GOP, GOP_len, idx_frame, psnr_list(idx_frame), ssim_list(idx_frame));

    for idx_B_frame = 1:GOP_len - 1
        idx_frame = decoding_map(idx_B_frame, 1);

        recons_block_non_skipped = MH_Residual_recons_IBI(y_cell_all(:, idx_frame), Phi_full, Psi, ...
            lambda, MH_win_radius, recons_seq(:, :, decoding_map(idx_B_frame, 2)), ...
            recons_seq(:, :, decoding_map(idx_B_frame, 3)));

        % MC, reconstruction for skipped blocks
        is_skipping_list = identifier_list_all(:, idx_frame);
        recons_block = recons_block_non_skipped;
        for idx_block = 1:block_num
            if (is_skipping_list(idx_block) == 0) % not skipped
                continue;
            elseif (is_skipping_list(idx_block) == 1) % ref to right frame
                recons_block(:, idx_block) = recons_block_seq ...
                    (:, idx_block, decoding_map(idx_B_frame, 3));
            elseif (is_skipping_list(idx_block) == 2) % ref to left frame
                recons_block(:, idx_block) = recons_block_seq ...
                    (:, idx_block, decoding_map(idx_B_frame, 2));
            else % ref to both left and right frames
                recons_block(:, idx_block) = ...
                    (recons_block_seq(:, idx_block, decoding_map(idx_B_frame, 2)) + ...
                    recons_block_seq(:, idx_block, decoding_map(idx_B_frame, 3))) / 2;
            end
        end

        recons_block(recons_block > 255) = 255;
        recons_block(recons_block < 0) = 0;
        recons_block_seq(:, :, idx_frame) = recons_block;
        recons_seq(:, :, idx_frame) = col2im(recons_block, [block_size, block_size], [height, width], 'distinct');

        % evaluate
        psnr_list(idx_frame) = Psnr(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
        ssim_list(idx_frame) = ssim(Y_seq(:, :, idx_frame), recons_seq(:, :, idx_frame));
        fprintf(' B Frame: %d, PSNR = %.2f dB, SSIM = %.4f \n', idx_frame, psnr_list(idx_frame), ssim_list(idx_frame));
    end
end
frame_num = idx_frame + 1;

%% 
fprintf('\n%s\n\tFrame Number: %d\n', file_name, frame_num);
fprintf('\tMean Skipping Rate: %.2f\n', mean(skipping_rate_list));
fprintf('\tMean Sampling Rate: %.4f\n', mean(sampling_rate_list(1:frame_num)));
fprintf('\tMean PSNR = %.2f dB \n', mean(psnr_list(1:frame_num)));
fprintf('\tMean SSIM = %.4f \n', mean(ssim_list(1:frame_num)));
