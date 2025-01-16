% SPG-L1 algorithm for convex optimization problem

% reference:
% E. van den Berg and M. P. Friedlander, “Probing the Pareto frontier for basis pursuit solutions,”
% SIAM J. Sci. Comput., vol. 31, no. 2, pp. 890–912, 2008, doi: 10.1137/080714488.

function recons_img_block = frameReconsSPGL1(y_cell, M_list, block_size, A_full)
% Input:
% y_cell: block measurements
% M_list: number of block measurements
% block_size: block size
% A_full: full size of measurement matrix or sensing matrix

% Output:
% recons_img_block: reconstruction result in block order


N = block_size^2; % block length in vector
block_num = length(M_list);

recons_img_block = zeros(N, block_num); % reconstruction result in column

opts = spgSetParms('verbosity', 0);
% reconstruct block by block
for idx_block = 1:block_num
    M = M_list(idx_block);
    if (M == 0) % empty block
        s = zeros(N, 1);
    elseif (M == N) % full measuring block
        s = A_full \ y_cell{idx_block};
    else % sparse block
        A = A_full(1:M, :);
        b = y_cell{idx_block};
        s = spg_bp(A, b, opts);
    end
    % result recording
    recons_img_block(:, idx_block) = s;
end

end
