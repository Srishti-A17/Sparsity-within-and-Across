function [X] = LIP_ellone(NablaF,param)
%
% solves the l1 minimization problem
% min_x f(x) + param.tau * || x ||_1,
% where f is as defined in the experiment in the letter
% 'Sparsity Within and Across Overlapping Groups', I. Bayram, 2017.
%
%%% input variables : 
% NablaF : handle for the gradient of f
% param.x : holds the estimate of X
% param.alpha : step parameter that should be less than 1/L
% param.tau : weight of the l_1 norm
% param.MAX_ITER : maximum number of iterations
%
%%% output variable : X
%
% Ilker Bayram, Istanbul Technical University, 2017,
% ibayram@itu.edu.tr

soft = @(x,z) max(abs(x) - z,0) .* (x ./ (abs(x) + 1e-10)); 

wb = waitbar(0,'l1 dereverberation');
for iter = 1:param.MAX_ITER,
    waitbar(iter/param.MAX_ITER,wb);
    
    % the Landweber step
    
    Z = param.x - param.alpha * NablaF(param.x);
    
    % threshold
    
    param.x = soft(Z,param.alpha * param.tau);
    
end
close(wb);

X = param.x;