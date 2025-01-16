function GOP_list = GOPListDecoder(GOP_len_init, tau_GOP_set, GOP_len_set, identifier_list_all)
% only the related params and identifiers are requireds

[block_num, frame_num] = size(identifier_list_all);
GOP_list = [GOP_len_init];

skipping_rate_list = sum(identifier_list_all ~= 0)' / block_num;

idx_frame = 1;
GOP_len = GOP_len_init;
while (idx_frame + GOP_len < frame_num)
    skipping_rate_GOP = mean(skipping_rate_list(idx_frame+1:idx_frame+GOP_len));
    idx_frame = idx_frame + GOP_len;

    if (skipping_rate_GOP < tau_GOP_set(1))
        GOP_len = GOP_len_set(1);
    elseif (skipping_rate_GOP < tau_GOP_set(2))
        GOP_len = GOP_len_set(2);
    elseif (skipping_rate_GOP < tau_GOP_set(3))
        GOP_len = GOP_len_set(3);
    else
        GOP_len = GOP_len_set(4);
    end
    GOP_list = [GOP_list; GOP_len];

end
end