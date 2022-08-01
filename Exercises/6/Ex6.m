clear 
close all hidden
dbstop if error
clc
%%
%%%%% PART 1 %%%%%
%% Q1
a = [1;-3;2].';% coefficients of (delayed) input x
b = [2;-1].';% coefficients of (delayed) output y
%
mySys = tf(a,b,-1,'variable','z^-1');
% -1: indicates unknown sampling period
% show transfer function using "mySys"
[poles,zeros] = pzmap(mySys);% plot pzmap with "pzmap(mySys)"
%
% stable because all poles inside unit circle; becomes obvious when looking
% at impulse response (IR): it shows a decay
%
% IIR filter because of recursion of output (denominator unequal 1)
%
%%% remark: ltiview is an old function. Use "linearSystemAnalyzer" instead
% ltiview('bode',mySys);% OLD
linearSystemAnalyzer({'pzmap';'impulse';'bode'},mySys);% new
%
% plot IR
% h = impulse(mySys,200);% calculate IR for 200 sampling periods; plot(h)
h = impulse(mySys);% number of periodsd determined by MATLAB automatically
%
figure;
N = [10;20;50;100];
[n,x,y,p] = deal(cell(numel(N),1));
% n: time axis
% x: input signal
% y: output signal
% p: handle for plot
for nn=1:numel(N)
    n{nn} = (0:1:N(nn));
    x{nn} = cos(2.*n{nn});
    y{nn} = conv(x{nn},h);
    p{nn} =subplot(4,1,nn); plot(y{nn});
end
linkaxes([p{1},p{2},p{3},p{4}],'xy')
%
% given a complex input signal: x[n] = A*z0^n, which is processed with a
% filter h, the output is given by the convolution between h and x, i.e.:
% y[n] = \sum_{k}(h[k]*x[n-k]) = \sum_{k}(h[k]*A*z0^{n-k}) = ...
% = A*z0^n \sum{k}(h[k]*z0^-k) = x[n]*H(z0), where H(z0) is the transfer
% function/ z transform of h evaluated at z0
% Hence, the output is similar to the input but the amplitude is modulated
% with H(z0);
% because: cos(2n) = (1/2)*exp(1i.*2n) + (1/2)*exp(-1i.*2n) there are
% actually two excitation signals with z0_{1} = exp(1i.*2) and z0_{2} = exp(-1i.*2)
%
% different input signal lengths lead to different responses because
% 1) already the input signal is always different: x[n]*H(z0) is always different 
% 2) conv(h,x) has always a different length: L_h + L_x -1
%% Q2
a = [1];
b = [1;2;1].';
mySys = tf(a,b,-1,'variable','z^-1');
[poles,zeros] = pzmap(mySys);
%%% remark: ltiview is an old function. Use "linearSystemAnalyzer" instead
% ltiview('bode',mySys);% OLD
%
% pole on unit-circle -> unstable
%
linearSystemAnalyzer({'pzmap';'impulse';'bode'},mySys);% new
h = impulse(mySys);
%% Q3
% stable FIR filter
alpha = 0.9;% "moving average" --> low pass filter
a = [1;alpha].';
b = [1];
mySys = tf(a,b,-1,'variable','z^-1');
[poles,zeros] = pzmap(mySys);
%%% remark: ltiview is an old function. Use "linearSystemAnalyzer" instead
% ltiview('bode',mySys);% OLD
linearSystemAnalyzer({'pzmap';'impulse';'bode'},mySys);% new
%
alpha = (-10:0.1:10).';
% since FIR filter and no poles alwyas stable, independent of alpha
A = numel(alpha);
counterUnstable = 0;
counterStable = 0;
for aa=1:A
    a = [1;alpha(aa)].';
    b = [1];
    mySys = tf(a,b,-1,'variable','z^-1');
    [poles,~] = pzmap(mySys);
    if abs(poles)>=1
        counterUnstable = counterUnstable + 1;
        unstableAlpha(counterUnstable) = alpha(aa);
        disp(['Filter becomes unstable at alpha:',num2str(unstableAlpha(counterUnstable))])
    else
        counterStable = counterStable + 1;
        stableAlpha(counterStable) = alpha(aa);
    end
    % disp(['Finished ',num2str(aa),'/',num2str(A)])
end
%
% setting alpha<0 --> highpass filter
%% Q4
% IIR filter that can become unstable
alpha = 0.9;
a = [1];
b = [1;alpha].';
mySys = tf(a,b,-1,'variable','z^-1');
[poles,zeros] = pzmap(mySys);
%%% remark: ltiview is an old function. Use "linearSystemAnalyzer" instead
% ltiview('bode',mySys);% OLD
linearSystemAnalyzer({'pzmap';'impulse';'bode'},mySys);% new
%
alpha = (-10:0.1:10).';
A = numel(alpha);
counterUnstable = 0;
counterStable = 0;
for aa=1:A
    a = [1];
    b = [1;alpha(aa)].';
    mySys = tf(a,b,-1,'variable','z^-1');
    [poles,~] = pzmap(mySys);
    if abs(poles)>=1
        counterUnstable = counterUnstable + 1;
        unstableAlpha(counterUnstable) = alpha(aa);
        disp(['Filter becomes unstable at alpha:',num2str(unstableAlpha(counterUnstable))])
    else
        counterStable = counterStable + 1;
        stableAlpha(counterStable) = alpha(aa);
    end
    % disp(['Finished ',num2str(aa),'/',num2str(A)])
end
%%
%%%%% PART 2 %%%%%
%% Q5
H1 = [2;3;4];
H2 = [3;4;5;6];
H = conv(H1,H2);
% convolving two polynomials is equivalent to multiplying two polynomials.
%% Q6
% remark: MATLAB: sinc(x) = sin(pi*x)/(pi*x)
% our convention: sinc(x) = sin(x)/x
% therefore, to use MATLAB's sinc function to get our results, write
% mySinc(x) = sinc(x/pi)
%
deltaT = 1;
fs = 1/deltaT;
n = (-50:deltaT:50).';
n = (-5000:deltaT:5000).';
temp = sinc((2.*n)./pi);
x = temp.^(2);
y = upsample(x,2);% K=2 --> compress frequency axis by factor = 2
fs_upsample = 2*fs;
z = exp(1i.*n).*x;% shift frequency spectrum by 1 unit to right 
%
mySpectrumAnalyzer(x,fs,'normalizedFreq');
mySpectrumAnalyzer(x,fs,'frequency');
mySpectrumAnalyzer(y,fs_upsample,'normalizedFreq');
mySpectrumAnalyzer(z,fs,'normalizedFreq');
%% Q7
%
% CTFT: continuous-time Fourier transform
% CTFT(cos(omega0*x)) = pi*delta(omega + omega0) + pi*delta(omega - omega0)
%
fs = 16000;
deltaT = 1/fs;
T = 1;
t = (0:deltaT:T-deltaT).';
f0 = 2e3;
f1 = 2*f0;
% minimum sampling freqeuncy >= 2*2f1
% maximum detectable frequency: fs
% spacing frquency axis: fs/nfft
%
x1 = cos(2*pi*f0.*t);
x2 = 5.*cos(2*pi*f1.*t);
x = x1 + x2;
%
mySpectrumAnalyzer(x,fs,'frequency');
%% Q8
x = [1;1;1;1];
w = [1;-1;1];
y_lin = conv(x,w);
y_cc = cconv(x,w,4);
%
% obtain a linear convolution output using cconv via zero padding
% pad both signals to length L_x + L_w -1
y_lin2 = cconv([x;0;0],[w;0;0;0],6);
%
w2 = [1;0;0;0];
y2_lin = conv(x,w2);
y2_cc = cconv(x,w2,4);
%
% there is no difference between these convolutions because not this
% corresponds to an operation as in overlap-save --> outlook for next
% exercise