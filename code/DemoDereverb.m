% Matlab code for the experiment in
% 'Sparsity Within and Across Overlapping Groups', I. Bayram, 2017

%% load the data

% variables to be loaded : 
% H : the impulse responses in the STFT domain
% X : STFT coefficients of the clean signal
% fs : sampling frequency
% win : window used in the STFT and ISTFT
% hop : hop-size for the STFT

load('../data/DereverbData'); 


%% produce observations in the STFT domain

% reverberant STFT coefficients
Y = upfirdn(X.',H.',1,1);
Y = Y.';

% add noise 
SNR = 5; % target SNR

n = randn(size(Y)) + 1i * randn(size(Y));
n = n * sqrt(numel(n)/ sum(abs(n(:)).^2));
sig = sqrt( sum(abs(Y(:)).^2)/(10^(SNR/10))/numel(Y) );

Y = Y + sig * n; % noisy observations

%% now set nablaf

Nf = @(x)NablaF(x,H.',Y.');

F = H.';
Y = Y.';
pf = upfirdn(F,conj(F(end:-1:1,:)),1,1);
L2 = max( max( abs( fft(pf,size(X,2) + size(F,1),1) )) ); % this is an approximation to L

ellone.alpha = 0.9/L2; % param.alpha <  1/L

%% reconstruction using the l1 norm
ellone.tau = 0.0112; % this specific value is found with a sweep search
ellone.MAX_ITER = 1000;
ellone.x = Y; % initialization

X1 = LIP_ellone(Nf,ellone);
X1 = X1.';
X1clip = X1(:,1:size(X,2));
snrl1 = snr(X(:),(X(:) - X1clip(:)))

%% reconstruction using the SWAG penalty

% form the weight matrix
S = 15; % size of each group
M = ones(S,S);
M(1:S+1:end) = 0;
NN = size(X,1) / S;
gam = 40;
SWAG.W = gam * kron(eye(NN),M); % the weight matrix
SWAG.alpha = ellone.alpha; 
SWAG.lambda = ellone.tau / 4; % non-optimal choice for lambda
SWAG.beta = SWAG.lambda * SWAG.alpha;
SWAG.eta = 1.9 / (1 + SWAG.beta * max(sum(SWAG.W)) ); % eta should satisfy the constraint in the manuscript
SWAG.iter = 5; % this is the number of inner iterations, namely K in Algorithm 1
SWAG.MAX_ITER = ellone.MAX_ITER;
SWAG.x = Y; % initialization

[XSWAG] = LIP_SWAG(Nf,SWAG);
XSWAG = XSWAG.';
XSWAGclip = XSWAG(:,1:size(X,2)); % adjust size for SNR computation
snrSWAG = snr(X(:),(X(:) - XSWAGclip(:)))


%% reconstruction using the modified penalty
% form the W matrix
S = size(X,1)/2; % the number of positive frequencies
ind = (0:S-1) * 0.75;
t = (ind.^4) .* exp(-ind);
t = 0.5 * t / sum(t);
gam = 900;
param.W = gam * toeplitz(t,t');% the weight matrix is concentrated around the diagonal and Toeplitz
param.alpha = ellone.alpha;
param.lambda = ellone.tau / 4; % same lambda value as SWAG
param.beta = param.lambda * param.alpha;
param.eta = 1.9 / (1 +  param.beta * max(sum(param.W)));
param.iter = SWAG.iter;
param.MAX_ITER = ellone.MAX_ITER;
param.x = Y; % initialization

[XProp] = LIP_SWAG(Nf,param);
XProp = XProp.';
XPropclip = XProp(:,1:size(X,2));
snrProp = snr(X(:),(X(:) - XPropclip(:)))

%% visualize the outputs

% STFT parameters for visualization

Fr = [0 3500]; % frequency range for STFT display
clim = [-60 0]; % amplitude range to display (in dB)

% Original signal
ttl = 'Original';
Norm = max(abs(X(:)));
figure;
DispSTFT(X/Norm, fs, length(win), hop, Fr, clim);
title(ttl);
print(strcat('../results/',ttl),'-dpng');

% Noisy and reverberant observation
ttl = 'Reverberant';
Y2 = Y.';
Y2 = Y2(:,1:size(X,2));
Norm = max(abs(Y2(:)));
figure;
DispSTFT(Y2/Norm, fs, length(win), hop, Fr, clim);
title(ttl);
print(strcat('../results/',ttl),'-dpng');

% Dereverbed with the l1 norm
ttl = 'Dereverbed-L1';
Norm = max(abs(X1(:)));
figure;
DispSTFT(X1clip/Norm, fs, length(win), hop, Fr, clim);
title(ttl);
print(strcat('../results/',ttl),'-dpng');

% Dereverbed with the SWAG penalty
ttl = 'Dereverbed-SWAG';
Norm = max(abs(XSWAG(:)));
figure;
DispSTFT(XSWAGclip/Norm, fs, length(win), hop, Fr, clim);
title(ttl);
print(strcat('../results/',ttl),'-dpng');

% Dereverbed with the proposed modification to the SWAG penalty
ttl = 'Dereverbed-Proposed';
Norm = max(abs(XProp(:)));
figure;
DispSTFT(XPropclip/Norm, fs, length(win), hop, Fr, clim);
title(ttl);
print(strcat('../results/',ttl),'-dpng');