function [Hx] = NablaF(x,F,Y)
% computes the gradient of the function
% 0.5 * || Y - F * x ||_2^2 at x
% where F is a band-wise convolution operator
%
% Ilker Bayram, 2017


S = size(Y,1);

Hx = upfirdn(x,F,1,1);
Hx = Hx(1:S,:) - Y;
Hx = upfirdn(Hx,conj(F(end:-1:1,:)),1,1);
Hx = Hx(end-S+1:end,:);