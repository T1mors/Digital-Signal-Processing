function y = myCircConv(u,h,M)
% u: input signal
% h: filter of length N
% M: DFT length
%
% y: output signal of length M x L
%
uBuffer = buffer(u,M,0);
L = size(uBuffer,2);
y = nan(M,L);
%
zerosToAppendToFilter = M-numel(h);
h_zp = fft([h;zeros(zerosToAppendToFilter,1)]);
for ll=1:L
    tempBuffer = uBuffer(:,ll);
    U = diag(fft(tempBuffer,M));
    %
    % since signal and filter are both real-valued, output is also real
    y(:,ll) = ifft(U*h_zp,M,'symmetric');
end
%
end