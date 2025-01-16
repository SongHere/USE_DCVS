% MH with Tikhonov Regularization

% reference
% C. Chen, E. W. Tramel, and J. E. Fowler,“Compressed-sensing recovery of images and video using multihypothesis predictions,”in 2011 Conf. Rec. 45th Asilomar Conf. Signals, Systems and Computers (ASILOMAR), Nov. 2011, pp. 1193–1198. doi: 10.1109/ACSSC.2011.6190204.

function w = MH_ABCS(y, Phi, H, lambda)

% MH
Y = repmat(y, 1, size(H, 2));
A = Phi * H;
Ga = diag(sum((Y - A).^2));

L = (A' * A + (lambda^2) * Ga);

if (rank(L) < size(L, 1))
    w = pinv(L) * (A' * y);
else
    w = L \ (A' * y);
end

end
