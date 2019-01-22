function [x] = mSWAG(z,param)
% function [x] = mSWAG(z,param)
%
% applies the projected gradient algorithm (PGA) 
% for a specified number of iterations to minimize
% 0.5 * | z - x |^2 + beta * ( |x|_1 + gamma * x' * W * x )
% 
% if z is a matrix, each column is treated as a different input
%
%%% input variables : 
% z : the variable to apply the prox. operator
% param.x : current estimate of x
% param.W : the weight matrix
% param.eta : the eta parameter in the manuscript
% param.beta : the beta parameter in the manuscript
% param.iter : number of iterations for PGA
%
% Ilker Bayram, Istanbul Technical University, 2017,
% ibayram@itu.edu.tr


if isfield(param,'x'),
    x = abs(param.x);
else
    x = abs(z);
end

az = abs(z);
sz = z ./ (az + 1e-10);

% the projected gradient algorithm
for iter = 1:param.iter,
    % the gradient step
    x = x - param.eta * ( (x - az) + param.beta * (param.W * x  + 1) );
    % the projection step
    x = max(x,0);
end

% correct the directions of x
x = x .* sz;