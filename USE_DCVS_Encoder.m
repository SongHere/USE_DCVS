% -----------------------------------------------------------------------
% The MATLAB programs of (encoder)
% Unified Sparsity Estimation-Based Distributed Compressive Video Sensing

% Author: Zhen Song
% The results within the paper are test on MATLAB R2018a
% -----------------------------------------------------------------------

%%
clear;
% clc;

addpath('.\data');
addpath('.\encoder');
%% data load
% [hall; container; foreman; coastguard; soccer]
% frame_num: 330, 300, 300, 300, 300
file_name = 'hall';
[Y_seq, u, v] = yuvRead([file_name, '_qcif.yuv'], 176, 144, 330);

% or
% [CW_Seq, news_cif_Y, football_ntsc_frames1-50_Y, FourPeople_1280x720_frames51-100_Y]
% file_name = 'CW_Seq';
% load([file_name, '.mat']);


Y_seq = double(Y_seq);
[height, width, frame_num] = size(Y_seq);

%% arguments configration
tau_En_skip_mean = 14; % threshold for block skipping

% GOP structure
% I B ... B P/I
% GOP size update
I_frame_GOP_period = 5;
GOP_len_set = [2, 4, 8, 16];
GOP_len_init = GOP_len_set(2);
tau_GOP_set = [0.1, 0.3, 0.6];

% _k: keyframe
% _nk: non-key frame

% Unified Sparsity Estimation
qf_k = 60; % quality factor
qf_nk = 30;

% Measurement allocation
alpha = 1.7; % measurement allocation factor

block_size = 8; % block size (only for 8*8)
M_init_k = 27; % initial number of measurements
M_init_nk = 18;
Mc_k = 32; % number of estimated coefficients
Mc_nk = 24;

% calculable params
N = block_size^2;
block_num = height * width / N;

%% matrix construction
% random measruement matrix
Phi_state = 0; % the random seed
randn('seed', Phi_state);
Phi_full = randn(N, N);
Phi_full = orth(Phi_full)';

% JPEG 8*8 luminance quantization table
q_table = JPEG88QTable();
[~, q_table, idx_col2zig] = zigzagScanning(q_table);
% q_table: JPEG 8*8 luminance quantization table in raster scanning order
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

idx_zig_k = idx_col2zig(2:Mc_k); % the first Mc DCT coefficients in the zig-zag scanning order
Lambda_k = I_N(idx_zig_k, :); % permutation matrix
Gamma_k = I_N .* Q_k; % threshold matrix
U_k = Lambda_k * (Gamma_k \ (Psi * P_k)); % USE matrix

idx_zig_nk = idx_col2zig(2:Mc_nk);
Lambda_nk = I_N(idx_zig_nk, :);
Gamma_nk = I_N .* Q_nk;
U_nk = Lambda_nk * (Gamma_nk \ (Psi * P_nk));

%% initialization
y_init_all = zeros(M_init_k, block_num, frame_num);

% data_stream: all measurements of non-skipped blocks
% (without the number of measurements of each block)
% identifier_list_all: 1 bit/block for P frame, 2 bits/block for B frame
% the related params: see USE_DCVE_Decoder.m
% are required to transmit to the decoder

data_stream = [];
sampling_rate_list = zeros(frame_num, 1);
identifier_list_all = zeros(block_num, frame_num);
skipping_rate_list = zeros(frame_num, 1);

time_encoding = zeros(frame_num, 1);

%% sampling
Phi_add_k = Phi_full(M_init_k+1:end, :);
Phi_add_nk = Phi_full(M_init_nk+1:end, :);

% adaptive sampling
% first frame, I frame
idx_frame = 1;

tic;
% BCS initial sampling
img_block = im2col(Y_seq(:, :, idx_frame), [block_size, block_size], 'distinct');
y_init = Phi_k * img_block;

% USE: Unified Sparsity Estimation
is_skipping_list = false(block_num, 1);
y_init_stream = y_init(:, ~is_skipping_list);
M_list = USE(y_init_stream, U_k, N, alpha);

% Supplementary sampling
y_stream_add = supplementarySampling(img_block(:, ~is_skipping_list), Phi_add_k, M_list-M_init_k);

time_encoding(idx_frame) = toc;

% transmitted data
data_stream = [data_stream; y_init_stream(:); y_stream_add];

% memory
% actually, only the initial measurements of the current GOP are required
y_init_all(:, :, idx_frame) = y_init;
GOP_len = GOP_len_init;

% infor
sampling_rate_list(idx_frame) = sum(M_list) / height / width;
fprintf('GOP: 1, Size = %d\n I Frame: %d, Skipping Rate = %.2f, Sampling Rate = %.4f; Mean ME time: %.4f ms\n', ...
    GOP_len, idx_frame, skipping_rate_list(idx_frame), ...
    sampling_rate_list(idx_frame), time_encoding(idx_frame)*1000);

idx_frame = 0;
GOP_num = 0;
while (idx_frame + GOP_len < frame_num)
    GOP_num = GOP_num + 1;

    tic;
    encoding_map = idxMidGet(1, GOP_len+1); % the GOP structure
    encoding_map = encoding_map + idx_frame;

    % BCS initial sampling
    idx_frame = encoding_map(1, 3);
    img_block = im2col(Y_seq(:, :, idx_frame), [block_size, block_size], 'distinct');
    y_init = Phi_k * img_block;

    % Skipping
    if (mod(GOP_num, I_frame_GOP_period) == 0) % I frame
        is_skipping_list = zeros(block_num, 1);
    else % P frame
        is_skipping_list = isSkipping(y_init-y_init_all(:, :, idx_frame-GOP_len), M_init_k, tau_En_skip_mean);
    end

    % USE
    y_init_stream = y_init(:, ~is_skipping_list);
    M_list = USE(y_init_stream, U_k, N, alpha);

    % Supplementary sampling
    y_stream_add = supplementarySampling(img_block(:, ~is_skipping_list), Phi_add_k, M_list-M_init_k);

    time_encoding(idx_frame) = toc;
    % transmitted data
    data_stream = [data_stream; y_init_stream(:); y_stream_add];
    identifier_list_all(:, idx_frame) = is_skipping_list;

    y_init_all(:, :, idx_frame) = y_init;

    % infor output
    skipping_rate_list(idx_frame) = sum(identifier_list_all(:, idx_frame) > 0) / block_num;
    sampling_rate_list(idx_frame) = sum(M_list) / height / width;
    fprintf('GOP: %d, Size = %d\n I/P Frame: %d, Skipping Rate = %.2f, Sampling Rate = %.4f; Encoding time: %.2f ms\n', ...
        GOP_num, GOP_len, idx_frame, skipping_rate_list(idx_frame), ...
        sampling_rate_list(idx_frame), time_encoding(idx_frame)*1000);

    % rest B frames
    for idx_B_frame = 1:GOP_len - 1
        tic;
        idx_frame = encoding_map(idx_B_frame, 1);
        img_block = im2col(Y_seq(:, :, idx_frame), [block_size, block_size], 'distinct');
        y_init = Phi_nk * img_block;

        % Skipping
        is_skipping_list_left = isSkipping(y_init-y_init_all(1:M_init_nk, :, encoding_map(idx_B_frame, 2)), ...
            M_init_nk, tau_En_skip_mean);
        is_skipping_list_right = isSkipping(y_init-y_init_all(1:M_init_nk, :, encoding_map(idx_B_frame, 3)), ...
            M_init_nk, tau_En_skip_mean);
        is_skipping_list = is_skipping_list_left * 2 + is_skipping_list_right;

        % USE
        y_init_stream = y_init(:, ~is_skipping_list);
        M_list = USE(y_init_stream, U_nk, N, alpha);

        % Supplementary sampling
        y_stream_add = supplementarySampling(img_block(:, ~is_skipping_list), Phi_add_nk, M_list-M_init_nk);
        
        time_encoding(idx_frame) = toc;
        
        % transmitted data
        data_stream = [data_stream; y_init_stream(:); y_stream_add];
        identifier_list_all(:, idx_frame) = is_skipping_list;

        y_init_all(1:M_init_nk, :, idx_frame) = y_init;

        % infor output
        skipping_rate_list(idx_frame) = sum(identifier_list_all(:, idx_frame) > 0) / block_num;
        sampling_rate_list(idx_frame) = sum(M_list) / height / width;
        fprintf(' B Frame: %d, Skipping Rate = %.2f, Sampling Rate = %.4f, Encoding time: %.2f ms\n', ...
            idx_frame, skipping_rate_list(idx_frame), sampling_rate_list(idx_frame), time_encoding(idx_frame)*1000);
    end

    skipping_rate_GOP = mean(skipping_rate_list(encoding_map(1, 2)+1:encoding_map(1, 3)));
    if (skipping_rate_GOP < tau_GOP_set(1))
        GOP_len = GOP_len_set(1);
    elseif (skipping_rate_GOP < tau_GOP_set(2))
        GOP_len = GOP_len_set(2);
    elseif (skipping_rate_GOP < tau_GOP_set(3))
        GOP_len = GOP_len_set(3);
    else
        GOP_len = GOP_len_set(4);
    end
end

frame_num = idx_frame + 1;
skipping_rate_list = skipping_rate_list(1:frame_num);
identifier_list_all = identifier_list_all(:, 1:frame_num);

%%
fprintf('\n%s \n Skipping Rate: %.2f, Sampling Rate: %.4f\n Encoding time: %.4f ms/frame\n', ...
    file_name, mean(skipping_rate_list), length(data_stream)/block_num/frame_num/N, mean(time_encoding)*1000);
