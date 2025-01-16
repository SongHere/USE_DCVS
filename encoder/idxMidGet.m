% function idx_en_order = gopEncodeOrderGetIBI(GOP_len)
%
% encoder_order = 1:GOP_len;
%
%
% end

function idx_mid_list = idxMidGet(idx_left, idx_right)

idx_mid_list = [];

if (idx_left + 1 == idx_right)
    return;
else
    idx_mid = (idx_right + idx_left) / 2;
    idx_mid_list = [idx_mid, idx_left, idx_right; ...
        idxMidGet(idx_left, idx_mid); idxMidGet(idx_mid, idx_right)];
end

end