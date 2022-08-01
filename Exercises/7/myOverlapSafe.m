function y = myOverlapSafe(u,h,M)
% u: input signal
% h: filter of length N
% M: DFT length = 2N
%
% y: output signal of length N x L
%
uBuffer = buffer(u,M,M/2);% buffered signal
L = size(uBuffer,2);% number of buffers
N = length(h);
h_zp = fft([h;zeros(M-N,1)],M);
%
y = nan(M/2,L);
K = [zeros(N,N),eye(N,N)];
%
for ll=1:L
    tempBlock = uBuffer(:,ll);
    U = diag(fft(tempBlock,M));
    y_dft = U*h_zp;
    y(:,ll) = K*ifft(y_dft,M,'symmetric');
end