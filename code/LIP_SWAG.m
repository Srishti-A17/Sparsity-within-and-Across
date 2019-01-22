function [X] = LIP_SWAG(NablaF,param)
%
% solves the minimization problem
% min_x f(x) + P_W(x),
% where f and P_W is as defined in the experiment in the letter
% 'Sparsity Within and Across Overlapping Groups', I. Bayram, 2017.
%
%%% input variables : 
% NablaF : handle for the gradient of f
% param.x : holds the estimate of X
% param.alpha : step parameter that should be less than 1/L
% param.lambda  : the lambda parameter in the letter
% param.beta  : beta parameter in the letter
% param.eta  : eta parameter in the letter
% param.iter :  number of inner iterations, namely K, in the letter
% param.MAX_ITER  : number of iterations
%
%%% output variable : X
%
% Ilker Bayram, Istanbul Technical University, 2017,
% ibayram@itu.edu.tr

soft = @(x,z) max(abs(x) - z,0) .* (x ./ (abs(x) + 1e-10)); 


wb = waitbar(0,'SWAG dereverberation');
for iter = 1:param.MAX_ITER,
    waitbar(iter/param.MAX_ITER,wb);
    % the Landweber step
    
    Z = param.x - param.alpha * NablaF(param.x);
    
    % descent iterations on the penalty function
    param.x = (param.x).';
    param.x = DenoiseSpSWAG(Z.',param);
    param.x = (param.x).';
end
close(wb);

X = param.x;