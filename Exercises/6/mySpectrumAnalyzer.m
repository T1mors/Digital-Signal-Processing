function [] = mySpectrumAnalyzer(x,varargin)
%      x: time-domain signal; column vector; size: L x 1
%     fs: sampling rate (in Hz); size: 1 x 1
%
if ~isempty(varargin)
    fs = varargin{1};
    variableString = varargin{2};
end
if 0
    nfft = 2^nextpow2(numel(x));
else
    nfft = numel(x);
end
X = fft(x,nfft);
% make correct frequency axis
mySpectrum = fftshift(X);
%
figure;
if strcmpi(variableString,'frequency')
    % deltaF = fs/nfft;
    freqAxis = (linspace(-fs/2,fs/2,nfft)).';
    plot(freqAxis,abs(mySpectrum))
elseif strcmpi(variableString,'normalizedFreq')
    % deltaOmega = (2*pi)/nfft;
    omegaAxis = (linspace(-pi,pi,nfft)).';
    plot(omegaAxis,abs(mySpectrum))
end
end