function [KL] = calculate_kl_divergence(p_weights,q_weights)
%CALCULATE_ENTROPY Summary of this function goes here
%   Detailed explanation goes here
    KL = sum(p_weights.*log2(p_weights/q_weights));
end
