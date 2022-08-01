%% main script to run Ex. 7
%
if 0
%% Task 1
% define signals
% nToConsiderForSignalU = n>=0&n<=12;
u_sig1 = (0.9).^(0:1:12).';
%
N = 12;
h_filt1 = ones(N,1);
%
M = 16;
% 
%% task 1.2a
y_myCConv = myCircConv(u_sig1,h_filt1,M);
%
%% task 1.2b
y_conv = conv(u_sig1,h_filt1);
sum(abs(y_myCConv(M-(M-N+1)+1:M)-y_conv(M-(M-N+1)+1:M)))% last M-N+1 samples
% but samples 9 to 16 are correct
%% task 1.2c
y_cconv = cconv(u_sig1,h_filt1,M);
%% task 1.2d
y_cconv24 = cconv(u_sig1,h_filt1,24);%M=13+12-1; numel(u):13, numel(h):12;
%
figure;
stem(y_conv); hold on
stem(y_cconv); hold on;
stem(y_cconv24); hold on;
%% task 1.2e
[u_sig2,fs] = audioread('speech.wav');
M = (16*fs)/1000;
h_filt2 = randn(M,1);
energyOfFilter = sum(abs(h_filt2).^(2));
h_filt2 = h_filt2./sqrt(energyOfFilter);% normalization
% sanity check: energyOfNormalizedFilter = sum(abs(h_filt2).^(2));
%
figure; periodogram(h_filt2, rectwin(length(h_filt2)),length(h_filt2),fs)
y_audio_cconv = myCircConv(u_sig2,h_filt2,M);
soundsc(y_audio_cconv(:),fs)
%
end
%% Task 2
%
%% Task 2.2
u_sig1 = (0.9).^(0:1:12).';
%
M = 32;
N = M/2;
h_filt1 = [ones(12,1);zeros(N-12,1)];
y_myOverlapSave = myOverlapSafe(u_sig1,h_filt1,M);%
%
%compare now
disp([y_conv(1:N),y_myOverlapSave,y_myCConv])
%
%% Task 2.3
y_audio2_myOverlapSave = myOverlapSafe(u_sig2,h_filt2,(32*fs/1000));%
soundsc(y_audio2_myOverlapSave(:),fs)