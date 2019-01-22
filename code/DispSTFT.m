function DispSTFT(c, fs, N, Hop, F, clim, varargin)
% Display the short-time Fourier transform coefficients c
%
% INPUT
%   c: STFT coefficients (2D array)
%   fs: sampling rate
%   N: length of FFT
%   Hop: Hop-size
%   F: display frequency range (in Hz)
%   clim: for 'Clim' of gca
%
% Ilker Bayram, Istanbul Technical University, 2017.

Fd = 1+round(F*N/fs);  

cdb = 20 * log10( abs( c(Fd(1):1:Fd(2),:) ) );

imagesc([0 size(c,2)*Hop/fs],F/1000,cdb);

set(gca,'Clim',clim);

axis xy;
xlabel( 'Time (seconds)' )
ylabel( 'Frequency (kHz)' )

colorbar