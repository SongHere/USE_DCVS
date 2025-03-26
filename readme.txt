USE_DCVS: 
 Unified Sparsity Estimation-based Distributed Compressive Video Sensing

The source codes and data of the work
Song, Z., Chen, J. Adaptive rate compression for distributed video sensing in wireless visual sensor networks. Vis Comput (2025). https://doi.org/10.1007/s00371-025-03845-5.

% -----------------------------------------------------------------------
% Author: Zhen Song
% The results within the paper are test on MATLAB R2018a
% -----------------------------------------------------------------------


Project structures:
 USE_DCVS_Encoder.m
 USE_DCVS_Decoder.m
 other dirs are related utilities


The decoder is implemented based on:
 the multi-hypothesis prediction[1] and residual reconstruction[2] techniques (reconstrution for non-key frames)
 the SPG-L1[3] algorithm (solving the corresponding convex optimization problems)

References
[1] C. Chen, E. W. Tramel, and J. E. Fowler, “Compressed-sensing recovery of images and video using multihypothesis predictions,” in 2011 Conference Record of the Forty Fifth Asilomar Conference on Signals, Systems and Computers (ASILOMAR), Nov. 2011, pp. 1193–1198. doi: 10.1109/ACSSC.2011.6190204.
[2] S. Mun and J. E. Fowler, “Residual Reconstruction for Block-Based Compressed Sensing of Video,” in 2011 Data Compression Conference, Mar. 2011, pp. 183–192. doi: 10.1109/DCC.2011.25.
[3] E. Van Den Berg and M. P. Friedlander, “Probing the Pareto Frontier for Basis Pursuit Solutions,” SIAM J. Sci. Comput., vol. 31, no. 2, pp. 890–912, Jan. 2009, doi: 10.1137/080714488.
https://my.ece.msstate.edu/faculty/fowler/BCSSPL/
https://friedlander.io/spgl1/
