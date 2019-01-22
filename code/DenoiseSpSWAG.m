function [X] = DenoiseSpSWAG(Y,param)
% function [x] = mSWAG(z,param)
%
% wrapper function for mSWAG, denoises the STFT
%
% applies the projected gradient algorithm (PGA) 
% for a specified number of iterations to minimize
% 0.5 * | z - x |^2 + beta * ( |x|_1 + gamma * x' * W * x )
% 
% if z is a matrix, each column is treated as a different input
%
%%% input variables : 
% z : the variable to threshold
% param.x : current estimate of x
% param.W : the weight matrix
% param.alpha : the alpha parameter
% param.lambda : the lambda parameter
% param.gamma : the gamma parameter
% param.beta : the step size
% param.iter : number of iterations for PGA
%
% Ilker Bayram, Istanbul Technical University, 2017,
% ibayram@itu.edu.tr


K = size(param.W,1);
Yrshp = reshape(Y, K, numel(Y) / K );
if isfield(param,'x'),
    param.x = reshape(param.x, K, numel(param.x) / K );
end
X = mSWAG(Yrshp, param);
X = reshape(X, size(Y,1), size(Y,2) );